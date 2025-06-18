"""
Read/write functionality for IGD.
"""

from contextlib import AbstractContextManager
from dataclasses import dataclass
from enum import Enum
from typing import BinaryIO, List, Tuple, Optional
import struct


class BpPosFlags(Enum):
    """
    An enumeration of bitwise flags that describe the contents of each variants genotype data.
    """

    MASK = 0xFF00000000000000
    SPARSE = 0x0100000000000000
    IS_MISSING = 0x0200000000000000


# Left shift amount to get the numCopies value
BP_POS_COPY_SHIFT = 48
# Mask for getting only the numCopies value
BP_POS_COPY_MASK = 0xFF << BP_POS_COPY_SHIFT
# Mask for getting only the bp position value
BP_POS_ONLY_MASK = ~(BpPosFlags.MASK.value | BP_POS_COPY_MASK)


def _read_uint64(file_obj: BinaryIO) -> int:
    return struct.unpack("Q", file_obj.read(8))[0]


def _read_uint32(file_obj: BinaryIO) -> int:
    return struct.unpack("I", file_obj.read(4))[0]


def _read_string(version: int, file_obj: BinaryIO) -> str:
    if version == 3:
        length = _read_uint64(file_obj)
    else:
        length = _read_uint32(file_obj)
    return file_obj.read(length).decode("utf-8")


def _read_u32_list(file_obj: BinaryIO) -> List[int]:
    result = []
    length = _read_uint32(file_obj)
    for _ in range(length):
        result.append(_read_uint32(file_obj))
    return result


def _div_round_up(value: int, round_to: int) -> int:
    return (value + (round_to - 1)) // round_to


def _samples_for_bv(data: bytes, index: int, sample_list: List[int]):
    mask = 0x1 << 7
    value = data[index]
    sample_offset = index * 8
    bit = 0
    while value != 0:
        is_set = 0 != (value & mask)
        if is_set:
            sample_list.append(sample_offset + bit)
        value = (value << 1) & 0xFF
        bit += 1


def _bv_for_bv(data: bytes, index: int, result_bv: List[int]):
    mask = 0x1 << 7
    value = data[index]
    sample_offset = index * 8
    bit = 0
    while value != 0:
        is_set = 0 != (value & mask)
        if is_set:
            result_bv[sample_offset + bit] = 1
        value = (value << 1) & 0xFF
        bit += 1


def flags_is_missing(flags: int):
    """
    Returns true if the flags specify that the variant represents missing data.

    :param flags: The flags, e.g. as returned from `get_position_and_flags`
    """
    return bool(flags & BpPosFlags.IS_MISSING.value)


# Internal constants shared between IGDReader and IGDWriter.
class IGDConstants:
    NUM_HEADER_BYTES = 128
    INDEX_ENTRY_BYTES = 16

    HEADER_FORMAT = "QQIIQIIQQQQQQQQQQQ"
    HEADER_MAGIC = 0x3A0C6FD7945A3481
    SUPPORTED_FILE_VERSION = 4

    FLAG_IS_PHASED = 0x1


class IGDReader:
    """
    Construct an IGDReader object for reading data from a file.

    :param file_obj: The file object to read from; should be opened in binary mode ("rb").
    """

    def __init__(self, file_obj: BinaryIO):
        self.file_obj = file_obj

        header_bytes = self.file_obj.read(IGDConstants.NUM_HEADER_BYTES)
        (
            magic,
            self._version,
            self._ploidy,
            _,
            self._num_var,
            self._num_idv,
            _,
            self._flags,
            self._fp_idx,
            self._fp_vars,
            self._fp_ind_ids,
            self._fp_var_ids,
            _,
            _,
            _,
            _,
            _,
            _,
        ) = struct.unpack(IGDConstants.HEADER_FORMAT, header_bytes)
        assert (
            magic == IGDConstants.HEADER_MAGIC
        ), "Invalid magic number; not an IGD file"
        assert (
            IGDConstants.SUPPORTED_FILE_VERSION == 4
        ), "When incrementing file version, check backwards compat below"
        assert self._version in (
            3,
            IGDConstants.SUPPORTED_FILE_VERSION,
        ), f"Unsupported IGD file format verison {self._version}"

        self._source = _read_string(self._version, self.file_obj)
        self._description = _read_string(self._version, self.file_obj)

        self._before_first_var = self.file_obj.tell()
        self._all_refs = None
        self._all_alts = None

    def _read_allele_info(self):
        self.file_obj.seek(self._fp_vars)
        self._all_refs = []
        self._all_alts = []
        for i in range(self.num_variants):
            ref = _read_string(self._version, self.file_obj)
            self._all_refs.append(ref)
            alt = _read_string(self._version, self.file_obj)
            self._all_alts.append(alt)

    @property
    def is_phased(self):
        """
        True if the data is phased. IGD doesn't support mixed phasedness.
        """
        return bool(self._flags & IGDConstants.FLAG_IS_PHASED)

    @property
    def version(self):
        """
        IGD file format version.
        """
        return self._version

    @property
    def ploidy(self):
        """
        Ploidy of each individual, between 1 and 8.
        """
        return self._ploidy

    @property
    def num_variants(self):
        """
        Number of variants in the file. This is not necessarily the same as the number of
        variants in the equivalent VCF file (for example), since IGD stores multi-allelic
        variants as multiple bi-allelic variants with the same base-pair position.
        """
        return self._num_var

    @property
    def num_individuals(self):
        """
        Number of individuals.
        """
        return self._num_idv

    @property
    def num_samples(self):
        """
        Number of samples. For phased data, this is num_individuals * ploidy. For unphased
        data this is just num_individuals. Every returned sample index (from get_samples())
        will be less than num_samples.
        """
        return self._num_idv * self._ploidy if self.is_phased else self._num_idv

    @property
    def source(self):
        """
        Source description of where the the IGD file came from.
        """
        return self._source

    @property
    def description(self):
        """
        Description of the IGD file.
        """
        return self._description

    def _get_var_idx_offset(self, variant_idx: int) -> int:
        return self._fp_idx + (variant_idx * IGDConstants.INDEX_ENTRY_BYTES)

    def get_position_flags_copies(self, variant_idx: int) -> Tuple[int, int, int]:
        """
        Given a variant index between 0...(num_variants-1), return the tuple (position, flags, num_copies).

        Much faster than `get_samples` or `get_samples_bv` because it only scans the variant index
        and does not read the actual genotype data.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The tuple (position, flags, num_copies) where position is the base-pair position
            (integer), flags is an integer that can be bitwise ANDed with BpPosFlags values, and
            num_copies is an integer indicating how many copies of the alternate allele this
            variant represents (for unphased data only).
        """
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        (bpp,) = struct.unpack(
            "Q", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES // 2)
        )
        flags = bpp & BpPosFlags.MASK.value
        copies = (bpp & BP_POS_COPY_MASK) >> BP_POS_COPY_SHIFT
        position = bpp & BP_POS_ONLY_MASK
        return (position, flags, copies)

    def get_position_and_flags(self, variant_idx: int) -> Tuple[int, int]:
        """
        Given a variant index between 0...(num_variants-1), return the tuple (position, flags).

        Much faster than `get_samples` or `get_samples_bv` because it only scans the variant index
        and does not read the actual genotype data.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The tuple (position, flags) where position is the base-pair position (integer) and
            flags is an integer that can be bitwise ANDed with BpPosFlags values.
        """
        assert self.is_phased, "Use get_position_flags_copies() for unphased data"
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        (bpp,) = struct.unpack(
            "Q", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES // 2)
        )
        flags = bpp & BpPosFlags.MASK.value
        position = bpp & BP_POS_ONLY_MASK
        return (position, flags)

    def get_samples(self, variant_idx: int) -> Tuple[int, bool, List[int]]:
        """
        Given a variant index between 0...(num_variants-1), return the tuple (position, missing, samples)
        where position is the base-pair position, missing is True when this represents a row of missing
        data, samples is a list of sample indexes that have the alt allele for the given variant, and
        copies is the number of copies of the alternate allele (for unphased data only).

        When missing is True, the sample list contains samples that are missing the given variant.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The tuple (position, is_missing, samples).
        """
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        bpp, fp_data = struct.unpack(
            "QQ", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES)
        )
        flags = bpp & BpPosFlags.MASK.value
        is_missing = 0 != (flags & BpPosFlags.IS_MISSING.value)
        is_sparse = 0 != (flags & BpPosFlags.SPARSE.value)
        position = bpp & BP_POS_ONLY_MASK
        self.file_obj.seek(fp_data)
        if is_sparse:
            sample_list = _read_u32_list(self.file_obj)
        else:
            sample_list = []
            byte_count = _div_round_up(self.num_samples, 8)
            data = self.file_obj.read(byte_count)
            for i in range(byte_count):
                _samples_for_bv(data, i, sample_list)
        return (position, is_missing, sample_list)

    def get_samples_bv(self, variant_idx: int):
        """
        Given a variant index between 0...(num_variants-1), return the tuple (position, missing, sample_bv)
        where position is the base-pair position, missing is True when this represents a row of missing
        data, and sample_bv is a bitvector representing samples that have the alt allele for the given variant.

        When missing is True, the sample vector contains samples that are missing the given variant.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The tuple (position, is_missing, samples).
        """
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        bpp, fp_data = struct.unpack(
            "QQ", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES)
        )
        flags = bpp & BpPosFlags.MASK.value
        is_missing = 0 != (flags & BpPosFlags.IS_MISSING.value)
        is_sparse = 0 != (flags & BpPosFlags.SPARSE.value)
        position = bpp & BP_POS_ONLY_MASK
        self.file_obj.seek(fp_data)
        byte_count = _div_round_up(self.num_samples, 8)
        sample_bv = [0 for _ in range(self.num_samples)]
        if is_sparse:
            as_list = _read_u32_list(self.file_obj)
            for sample_idx in as_list:
                sample_bv[sample_idx] = 1
        else:
            data = self.file_obj.read(byte_count)
            for i in range(byte_count):
                _bv_for_bv(data, i, sample_bv)
        return (position, is_missing, sample_bv)

    def get_alt_allele(self, variant_idx: int) -> str:
        """
        Get the alternative allele string for the given variant index.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The string representation of the alternate allele.
        """
        if self._all_alts is None:
            self._read_allele_info()
            assert self._all_alts is not None  # Make mypy happy
        return self._all_alts[variant_idx]

    def get_ref_allele(self, variant_idx: int) -> str:
        """
        Get the reference allele string for the given variant index.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The string representation of the reference allele.
        """
        if self._all_refs is None:
            self._read_allele_info()
            assert self._all_refs is not None  # Make mypy happy
        return self._all_refs[variant_idx]

    def get_individual_ids(self) -> List[str]:
        """
        Get a list of identifiers for the individuals in this dataset. The 0th individual's label
        is at list position 0, and the last individual is at list position (num_individuals-1).

        :return: Empty list if there are no identifiers (it is optional). Otherwise a list of strings.
        """
        result = []
        if self._fp_ind_ids > 0:
            self.file_obj.seek(self._fp_ind_ids)
            count = _read_uint64(self.file_obj)
            assert (
                count == self._num_idv
            ), f"Malformed file: {count} individual labels but expected {self._num_idv}"
            for _ in range(count):
                result.append(_read_string(self._version, self.file_obj))
        return result

    @property
    def has_variant_ids(self) -> bool:
        """
        True if there are variant IDs in this IGD file.
        """
        return self._fp_var_ids > 0

    def get_variant_ids(self) -> List[str]:
        """
        Get a list of identifiers for the variants in this dataset. The 0th variants's label
        is at list position 0, and the last variant is at list position (num_variants-1).

        :return: Empty list if there are no identifiers (it is optional). Otherwise a list of strings.
        """
        result = []
        if self._fp_var_ids > 0:
            self.file_obj.seek(self._fp_var_ids)
            count = _read_uint64(self.file_obj)
            assert (
                count == self._num_var
            ), f"Malformed file: {count} variant labels but expected {self._num_var}"
            for _ in range(count):
                result.append(_read_string(self._version, self.file_obj))
        return result

    def lower_bound_position(self, position) -> int:
        """
        Return the first variant index with position that is greater than or equal to the given position.
        Will return `num_variants` if the given position is greater than all positions in the IGD.

        :param position: The position to search for.
        :type position: int
        :return: The first variant index with position greater than or equal to the given position.
        :rtype: int
        """
        low = 0
        high = self.num_variants - 1
        mid = high
        while low <= high:
            mid = low + ((high - low) // 2)
            mid_pos, _ = self.get_position_and_flags(mid)
            if mid_pos < position:
                low = mid + 1
            elif mid_pos > position:
                high = mid - 1
            else:
                return mid
        return low


class IGDFile(IGDReader, AbstractContextManager):
    """
    DEPRECATED DO NOT USE. Context object for loading an IGD file. See IGDReader instead.

    :param filename: The filename to open.
    """

    def __init__(self, filename: str):
        file_obj = open(filename, "rb")
        super().__init__(file_obj)

    def __del__(self):
        if self.file_obj is not None:
            self.file_obj.close()
        self.file_obj = None

    def __exit__(self, *args):
        if self.file_obj is not None:
            self.file_obj.close()
        self.file_obj = None


# Internal class for managing the fixed-sized header of an IGD file.
@dataclass
class IGDHeader:
    magic: int
    version: int
    ploidy: int
    sparse_threshold: int
    num_variants: int
    num_individuals: int
    reserved: int
    flags: int
    fp_index: int
    fp_variants: int
    fp_individualids: int
    fp_variantids: int

    def pack(self) -> bytes:
        return struct.pack(
            IGDConstants.HEADER_FORMAT,
            self.magic,
            self.version,
            self.ploidy,
            self.sparse_threshold,
            self.num_variants,
            self.num_individuals,
            0,
            self.flags,
            self.fp_index,
            self.fp_variants,
            self.fp_individualids,
            self.fp_variantids,
            0,
            0,
            0,
            0,
            0,
            0,
        )


def _write_u64(file_obj, value):
    file_obj.write(struct.pack("Q", value))


def _write_u32(file_obj, value):
    file_obj.write(struct.pack("I", value))


def _write_string(file_obj, string):
    _write_u32(file_obj, len(string))
    file_obj.write(string.encode("utf-8"))


class IGDWriter:
    """
    Construct an IGDWriter for a given output stream.

    :param out_stream: The output stream to write to; usually a file opened via mode "wb".
    :param individuals: The number of individual samples in the file. NOT the number of haploids
        unless ploidy=1.
    :param ploidy: The ploidy of each individual sample.
    :param phased: Whether the data being stored is phased.
    :param source: A string describing where the data came from.
    :param description: A string describing the contents of the file.
    :param sparse_threshold: The threshold for choosing between sparse and dense sample lists when
        writing variant data to the file. Default is 32, which means that we still store variants
        sparsely if their frequency is less than or equal to 1/32. This is the threshold that is
        theoretically the break-even point between sparse and dense representations (since the
        sparse representation uses 32-bit integers, and dense uses a bit per sample).
    """

    def __init__(
        self,
        out_stream: BinaryIO,
        individuals: int,
        ploidy: int = 2,
        phased: bool = True,
        source: str = "",
        description: str = "",
        sparse_threshold: int = 32,
    ):
        self.out = out_stream
        self.header = IGDHeader(
            magic=IGDConstants.HEADER_MAGIC,
            version=IGDConstants.SUPPORTED_FILE_VERSION,
            ploidy=ploidy,
            sparse_threshold=sparse_threshold,
            num_variants=0,
            num_individuals=individuals,
            reserved=0,
            flags=0 if not phased else IGDConstants.FLAG_IS_PHASED,
            fp_index=0,
            fp_variants=0,
            fp_individualids=0,
            fp_variantids=0,
        )
        self.source = source
        self.description = description
        self.ref_alleles: List[str] = []
        self.alt_alleles: List[str] = []
        self.index: List[bytes] = []
        self.num_samples = (individuals * ploidy) if phased else individuals
        self.should_be_sparse = self.num_samples / sparse_threshold
        # The position of the last variant that we wrote. These have to be in order.
        self.last_var_position = 0

    def write_header(self):
        """
        Write the file header to the current output buffer position. Fails if that
        position if not the start of the buffer.
        """
        assert self.out.tell() == 0, "Writing header to wrong location"
        self.out.write(self.header.pack())
        _write_string(self.out, self.source)
        _write_string(self.out, self.description)

    @staticmethod
    def _make_index_entry(
        position: int, is_missing: bool, is_sparse: bool, num_copies: int, filepos: int
    ) -> bytes:
        flags = 0
        if is_sparse:
            flags |= BpPosFlags.SPARSE.value
        if is_missing:
            flags |= BpPosFlags.IS_MISSING.value
        encoded_pos = position | flags | (num_copies << BP_POS_COPY_SHIFT)
        return struct.pack("QQ", encoded_pos, filepos)

    def write_variant(
        self,
        position: int,
        ref_allele: str,
        alt_allele: str,
        samples: List[int],
        is_missing: bool = False,
        num_copies: int = 0,
    ):
        """
        Write the next variant, including sample information, to the file.
        Variants must be written in ascending order of their base pair position.

        :param position: Base-pair position.
        :param ref_allele: The reference allele.
        :param alt_allele: The alternate allele.
        :param samples: The list of samples, as indexes. E.g. the list [4, 10] means that samples
            numbered 4 and 10 have this variant's alternate allele. This list must be in ascending
            order.
        :param is_missing: [Optional] Set to true if the sample list represents the list of samples
            that are missing allele values at this position (in which case the reference and alt
            allele are somewhat irrelevant).
        :param num_copies: [Optional] For unphased data, set to the number of copies of the alternate
            allele that this variant represents (1...ploidy).
        """
        assert (
            position >= self.last_var_position
        ), "Out of order variant (must be written in ascending order)"
        is_phased = self.header.flags & IGDConstants.FLAG_IS_PHASED
        assert (
            is_phased or num_copies > 0
        ), "Unphased data must be written with num_copies specified to a non-zero value"
        assert num_copies == self.header.ploidy or not (
            not is_phased and is_missing
        ), "Unphased missing rows must have num_copies equal to ploidy (missingness is for an entire individual)"
        self.last_var_position = position
        self.ref_alleles.append(ref_allele)
        self.alt_alleles.append(alt_allele)
        is_sparse = len(samples) <= self.should_be_sparse
        filepos = self.out.tell()
        self.index.append(
            self._make_index_entry(position, is_missing, is_sparse, num_copies, filepos)
        )
        if is_sparse:
            _write_u32(self.out, len(samples))
            prev_sample = -1
            for sample_idx in samples:
                assert (
                    prev_sample < sample_idx
                ), "Sample indexes must be in ascending order"
                assert sample_idx < self.num_samples, "Invalid sample index"
                _write_u32(self.out, sample_idx)
                prev_sample = sample_idx
        else:
            num_bytes = (self.num_samples + 7) // 8
            # This is not particularly efficient, as we make two passes over the data.
            # If you want max efficiency, you should probably be using the C++ version.
            data = [0 for _ in range(num_bytes)]
            for sample_idx in samples:
                assert sample_idx < self.num_samples, "Invalid sample index"
                element = sample_idx // 8
                bit = 7 - (sample_idx % 8)
                data[element] = data[element] | (1 << bit)
            self.out.write(struct.pack("B" * num_bytes, *data))
        self.header.num_variants += 1

    def write_index(self):
        """
        Write the variant position index. Must be called _after_  all calls to write_variant().
        """
        assert len(self.index) == self.header.num_variants
        self.header.fp_index = self.out.tell()
        for item in self.index:
            self.out.write(item)

    def write_variant_info(self):
        """
        Write the variant information table. Must be called _after_ all calls to write_variant().
        """
        assert len(self.ref_alleles) == self.header.num_variants
        assert len(self.alt_alleles) == self.header.num_variants
        self.header.fp_variants = self.out.tell()
        for i in range(self.header.num_variants):
            _write_string(self.out, self.ref_alleles[i])
            _write_string(self.out, self.alt_alleles[i])

    def write_individual_ids(self, labels: List[str]):
        """
        Write the identifiers for the individual samples (optional).

        :param labels: Empty list or a list of strings, one for each individual.
        """
        if len(labels) == 0:
            self.header.fp_individualids = 0
        else:
            assert len(labels) == self.header.num_individuals
            self.header.fp_individualids = self.out.tell()
            _write_u64(self.out, len(labels))
            for label in labels:
                _write_string(self.out, label)

    def write_variant_ids(self, labels: List[str]):
        """
        Write the identifiers for the variants (optional).

        :param labels: Empty list or a list of strings, one for each variant. That is, if you called
            write_variant() `X` times, then there should be `X` entries in this list.
        """
        if len(labels) == 0:
            self.header.fp_variantids = 0
        else:
            assert len(labels) == self.header.num_variants
            self.header.fp_variantids = self.out.tell()
            _write_u64(self.out, len(labels))
            for label in labels:
                _write_string(self.out, label)


class IGDTransformer:
    """
    Class for transforming one IGD file to another.

    :param in_stream: The input stream for the input IGD file. Usually a file opened via
        mode "rb".
    :param out_stream: The output stream for the output IGD file. Usually a file opened via
        mode "wb".
    :param use_bitvectors: If True, the modify_samples callback will be invoked with a
        a List if 1s and 0s, where position "i" being 1 means sample "i" has the alternate allele.
    """

    def __init__(
        self, in_stream: BinaryIO, out_stream: BinaryIO, use_bitvectors: bool = False
    ):
        self.reader = IGDReader(in_stream)
        self.writer = IGDWriter(
            out_stream,
            self.reader.num_individuals,
            self.reader.ploidy,
            self.reader.is_phased,
            self.reader.source,
            self.reader.description,
        )
        self.use_bitvectors = use_bitvectors

    def transform(self):
        """
        Transform the input file to the output file, invoking modify_samples() for every variant
        from the input file. If modify_samples() returns None then the variant will not be emitted
        to the output file. Otherwise the variant will be emitted with whatever sample list is
        returned from modify_samples().
        """
        variant_ids = self.reader.get_variant_ids()
        self.writer.write_header()
        for i in range(self.reader.num_variants):
            _, _, num_copies = self.reader.get_position_flags_copies(i)
            if self.use_bitvectors:
                position, is_missing, samples = self.reader.get_samples_bv(i)
            else:
                position, is_missing, samples = self.reader.get_samples(i)
            samples = self.modify_samples(position, is_missing, samples, num_copies)
            if samples is None:
                variant_ids[i] = None
            else:
                # If the user is using bitvectors convert back to a sample list before writing.
                if self.use_bitvectors:
                    samples_list = []
                    for s in range(len(samples)):
                        if samples[s]:
                            samples_list.append(s)
                    samples = samples_list
                self.writer.write_variant(
                    position,
                    self.reader.get_ref_allele(i),
                    self.reader.get_alt_allele(i),
                    samples,
                    is_missing,
                    num_copies,
                )
        self.writer.write_index()
        self.writer.write_variant_info()
        self.writer.write_individual_ids(self.reader.get_individual_ids())
        self.writer.write_variant_ids(
            list(filter(lambda v: v is not None, variant_ids))
        )
        self.writer.out.seek(0)
        self.writer.write_header()

    def modify_samples(
        self,
        position: int,
        is_missing: bool,
        samples: List[int],
        num_copies: int = 0,
    ) -> Optional[List[int]]:
        raise NotImplementedError("Derived class must implement modify_samples()")
