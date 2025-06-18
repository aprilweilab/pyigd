"""
Extra functionality that is not core to accessing an IGD, but helpful for manipulating
the information in/related to an IGD.
"""

from pyigd.readwrite import (
    _div_round_up,
    _read_uint32,
    _read_string,
    _write_string,
    _write_u64,
    BpPosFlags,
    IGDConstants,
    IGDReader,
    IGDWriter,
    flags_is_missing,
)
from typing import List, BinaryIO, Tuple, Optional
import os
import struct

try:
    import numpy
except ImportError:
    numpy = None  # type: ignore


def collect_next_site(igd_reader, variant_index: int) -> List[int]:
    """
    Given a variant index to start at, iterate all consecutive variants that have the same position
    and return the variant indices for that position. The returned indices are in ascending order.

    :param igd_reader: The IGDReader representing the IGD file.
    :type igd_reader: IGDReader
    :param variant_index: The variant index to start scanning from.
    :type variant_index: int
    :return: The list of subsequent variant indexes that all share the same position.
    :rtype: List[int]
    """
    results = []
    site_position = None
    while variant_index < igd_reader.num_variants:
        position, _, _ = igd_reader.get_position_flags_copies(variant_index)
        if site_position is not None and position != site_position:
            break
        site_position = position
        results.append(variant_index)
        variant_index += 1
    return results


# For unphased data, we can have an individual which has 1 copy of A, 2 copies of G, and a ploidy of 4.
# Then the reference has 1 copy (the inverse). We need to look @ all copies to determine ploidy copies
def get_inverse_sample_list(igd_reader, variant_indices: List[int]) -> List[List[int]]:
    """
    Given a list of variant indices, compute the set of samples that are covered by those indices, and then
    invert that list. Works for both phased and unphased data. One usage for this function is to get an explicit
    list of samples that have the reference, since the reference is stored implicitly in an IGD.

    If the input data has overlapping variants that add up to more than PLOIDY for some samples, then
    we return an empty list. In real data this can happen at sites containing indels/SVs, depending
    on how the dataset was created. If you don't want this behavior (skipping the whole site) then
    you must filter out these problematic variants in the IGD file.

    :param igd_reader: The IGDReader representing the IGD file.
    :type igd_reader: IGDReader
    :param variant_indices: The indices of the variants who total sample set (i.e., all of them unioned together)
        you want the inverse of.
    :type variant_indices: List[int]
    :return: A list of sample lists, one for each number of copies between 1...PLOIDY. For phased data, there
        will always be a single sample list. For unphased data, there will always be PLOIDY sample lists.
    :rtype: List[List[int]]
    """
    # If you are inverting a sample list, then sparseness no longer matters. Either:
    # - The original list was not sparse, so you have to iterate a lot of items.
    # - The original list was sparse, so the inverse is not.
    # Hence, the inversion is a non-sparse operation and we use a matrix here.
    assert numpy is not None, "numpy required for this operation: 'pip install numpy'"
    ploidy = igd_reader.ploidy
    coverage_count = numpy.full(
        igd_reader.num_samples, 1 if igd_reader.is_phased else ploidy, dtype=numpy.int32
    )
    for i in variant_indices:
        _, flags, num_copies = igd_reader.get_position_flags_copies(i)
        if num_copies == 0:
            # For phased data.
            if igd_reader.is_phased:
                num_copies = 1
            # We could also assert this, but more robust to just force it to happen.
            elif flags_is_missing(flags):
                num_copies = ploidy
        _, _, samples = igd_reader.get_samples(i)
        this_vector = numpy.zeros(igd_reader.num_samples, dtype=numpy.int32)
        numpy.put(this_vector, samples, num_copies)
        coverage_count -= this_vector
    result = []
    # If we have negative values, it means that the input data had overlapping variants that added up to
    # more than PLOIDY for some samples. In real data this often happens at sites containing indels/SVs,
    # where the overlap isn't very well resolved and so a particular sample is called as having SNV "x" and
    # indel "y", which isn't always very interpretable.
    # We return no results for any such sites.
    if numpy.all(coverage_count >= 0):
        if igd_reader.is_phased:
            result.append(numpy.flatnonzero(coverage_count).tolist())
        else:
            for num_copies in range(1, ploidy + 1):
                result.append(numpy.flatnonzero(coverage_count == num_copies).tolist())
    return result


def igd_merge(
    out_file: str,
    in_readers: List[IGDReader],
    force_overwrite: bool = False,
    description: Optional[str] = None,
):
    """
    Merge the IGD files given by a list of IGDReader into a single output file with the given
    filename. The input files must be _mutually exclusive_ by genome range, such that if one
    input covers variants over the range (R1, R2) and another covers (R3, R4) then either
    R1 >= R4 or R3 >= R2. The samples described by the inputs must also be identical.

    If only some input IGDs have variant IDs, the remaining variants will get the empty string
    for their identifiers.
    The individual IDs will be used from the first (ascending order genetic position) input IGD
    and will not be checked against other IGD files individual IDs.

    :param out_file: The filename to write the output IGD to.
    :type out_file: str
    :param in_readers: List of IGDReader objects for the input IGDs to be merged.
    :type in_readers: List[IGDReader]
    :param force_overwrite: Optional. Set to True to always write the output file, even if it
        already exists.
    :type force_overwrite: bool
    :param description: Optional. Description to write to the IGD header. If not specified (None)
        then the description of the first input IGD will be used.
    :type description: Optional[str]
    """
    assert (
        not os.path.exists(out_file) or force_overwrite
    ), "Output file already exists. Use force_overwrite=True or remove the file."
    assert len(in_readers) >= 2, "No point in merging fewer than 2 IGD files"

    # Sanity check that properties of inputs are consistent.
    ploidy = in_readers[0].ploidy
    num_indiv = in_readers[0].num_individuals
    phased = in_readers[0].is_phased
    source = f"igd_merge({in_readers[0].source})"
    description = in_readers[0].description if description is None else description
    write_variant_ids = False
    for reader in in_readers:
        assert reader.ploidy == ploidy, "Multiple ploidy values in input IGDs"
        assert (
            reader.num_individuals == num_indiv
        ), "Different sample set sizes in input IGDs"
        assert reader.is_phased == phased, "Different phasedness in input IGDs"
        assert reader.num_variants > 0, "Empty IGD provided (no variants)"
        assert reader.version > 3, "Not support for old IGD file format versions"
        if reader.has_variant_ids:
            write_variant_ids = True

    # Now sort the input readers according to their first variant, and then ensure no overlap.
    in_readers.sort(key=lambda r: r.get_position_flags_copies(0)[0])
    for i in range(1, len(in_readers)):
        pr = in_readers[i - 1]
        r = in_readers[i]
        prev_end = pr.get_position_flags_copies(pr.num_variants - 1)
        curr_start = r.get_position_flags_copies(0)
        assert (
            prev_end <= curr_start
        ), f"Overlapping input IGDs: one ends at {prev_end} and the other starts at {curr_start}"

    def get_offset_for_variant(reader: IGDReader, index: int) -> Tuple[int, int]:
        end_offset = False
        if index == reader.num_variants:
            index = reader.num_variants - 1
            end_offset = True
        reader.file_obj.seek(reader._get_var_idx_offset(index))
        bpp, start_offset = struct.unpack(
            "QQ", reader.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES)
        )
        if end_offset:
            flags = bpp & BpPosFlags.MASK.value
            is_sparse = 0 != (flags & BpPosFlags.SPARSE.value)
            if is_sparse:
                byte_count = _read_uint32(reader.file_obj) * 4
            else:
                byte_count = _div_round_up(reader.num_samples, 8)
            return start_offset + byte_count, bpp
        return start_offset, bpp

    def copy_bytes(
        in_stream: BinaryIO,
        out_stream: BinaryIO,
        byte_size: int,
        max_chunk: int = 64 * 1024 * 1024,
    ):
        written = 0
        while written < byte_size:
            data = in_stream.read(max_chunk)
            out_stream.write(data)
            written += len(data)

    with open(out_file, "wb") as fout:
        writer = IGDWriter(fout, num_indiv, ploidy, phased, source, description)

        # We interpret as little as possible from the input IGDs, to speed up the merging.
        # The order of the resulting file looks like this:
        #   HEADER
        writer.write_header()

        #   SAMPLE LISTS (input 1)   <-- offset1
        #   ...
        #   SAMPLE LISTS (input N)   <-- offsetN
        old_offsets = []
        new_offsets = []
        for r in in_readers:
            start_offset, _ = get_offset_for_variant(r, 0)
            end_offset, _ = get_offset_for_variant(r, r.num_variants)
            assert end_offset > start_offset
            byte_size = end_offset - start_offset
            r.file_obj.seek(start_offset)
            old_offsets.append(start_offset)
            new_offsets.append(writer.out.tell())
            copy_bytes(r.file_obj, writer.out, byte_size)
            writer.header.num_variants += r.num_variants

        #   INDEX + offset1 (input1)
        #   ...
        #   INDEX + offsetN (inputN)
        writer.header.fp_index = writer.out.tell()
        for r, old, new in zip(in_readers, old_offsets, new_offsets):
            for i in range(r.num_variants):
                reader_offset, encoded_position = get_offset_for_variant(r, i)
                writer.out.write(
                    struct.pack("QQ", encoded_position, (reader_offset - old) + new)
                )

        #   VARIANT INFO (input1)
        #   ...
        #   VARIANT INFO (inputN)
        writer.header.fp_variants = writer.out.tell()
        for r in in_readers:
            r.file_obj.seek(r._fp_vars)
            for i in range(r.num_variants):
                ref = _read_string(r._version, r.file_obj)  # REF
                alt = _read_string(r._version, r.file_obj)  # ALT
                _write_string(writer.out, ref)
                _write_string(writer.out, alt)

        #   INDIVIDUAL IDS (input1) -- optional
        indiv_ids = in_readers[0].get_individual_ids()
        writer.write_individual_ids(indiv_ids)

        #   VARIANT IDS (input1)
        #   ...
        #   VARIANT IDS (input2)
        if write_variant_ids:
            writer.header.fp_variantids = writer.out.tell()
            _write_u64(writer.out, writer.header.num_variants)
            for r in in_readers:
                vids = r.get_variant_ids()
                if len(vids) < r.num_variants:
                    vids.extend([""] * (r.num_variants - len(vids)))
                for label in vids:
                    _write_string(writer.out, label)

        writer.out.seek(0)
        writer.write_header()
