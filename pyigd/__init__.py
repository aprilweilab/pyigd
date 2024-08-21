"""
Module for reading Indexable Genotype Data (IGD) files.
"""
import struct
from enum import Enum
from typing import BinaryIO, List, Tuple, Union
from dataclasses import dataclass
# Optional import. BitVector is only needed when using the APIs that return them.
try:
    import BitVector  # type: ignore
    BVType = BitVector.BitVector
except ImportError:
    BitVector = None
    BVType = List[int]


class BpPosFlags(Enum):
    """
    An enumeration of bitwise flags that describe the contents of each variants genotype data.
    """
    MASK = 0xFF00000000000000
    SPARSE = 0x0100000000000000
    IS_MISSING = 0x0200000000000000


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


def samples_for_bv(data: bytes, index: int, sample_list: List[int]):
    mask = 0x1 << 7
    value = data[index]
    sample_offset = index * 8
    bit = 0
    while value != 0:
        is_set = (0 != (value & mask))
        if is_set:
            sample_list.append(sample_offset + bit)
        value = (value << 1) & 0xFF
        bit += 1


class IGDConstants:
    NUM_HEADER_BYTES = 128
    INDEX_ENTRY_BYTES = 16

    HEADER_FORMAT = "QQIIQQQQQQQQQQQQQ"
    HEADER_MAGIC = 0x3a0c6fd7945a3481
    SUPPORTED_FILE_VERSION = 4

    FLAG_IS_PHASED = 0x1


class IGDReader:
    def __init__(self, file_obj: BinaryIO):
        """
        Construct an IGDReader object for reading data from a file.

        :param filename: The path to the file to be opened.
        """
        self.file_obj = file_obj

        header_bytes = self.file_obj.read(IGDConstants.NUM_HEADER_BYTES)
        (magic, self._version, self._ploidy, _, self._num_var, self._num_idv,
         self._flags, self._fp_idx, self._fp_vars, self._fp_ind_ids, self._fp_var_ids, _, _, _, _,
         _, _) = struct.unpack(IGDConstants.HEADER_FORMAT, header_bytes)
        assert magic == IGDConstants.HEADER_MAGIC, "Invalid magic number; not an IGD file"
        assert IGDConstants.SUPPORTED_FILE_VERSION == 4, "When incrementing file version, check backwards compat below"
        assert self._version in (3, IGDConstants.SUPPORTED_FILE_VERSION), f"Unsupported IGD file format verison {self._version}"

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
        Number of haploid samples (same as num_individuals * ploidy).
        """
        return self._num_idv * self._ploidy

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
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        bpp, _ = struct.unpack("QQ", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES))
        flags = (bpp & BpPosFlags.MASK.value)
        position = bpp & ~BpPosFlags.MASK.value
        return (position, flags)

    def get_samples(self, variant_idx: int) -> Tuple[int, bool, List[int]]:
        """
        Given a variant index between 0...(num_variants-1), return the tuple (position, missing, samples)
        where position is the base-pair position, missing is True when this represents a row of missing
        data, and samples is a list of sample indexes that have the alt allele for the given variant.

        When missing is True, the sample list contains samples that are missing the given variant.

        :param variant_idx: Variant index between 0...(num_variants-1). Variants are ordered from
            smallest to largest base-pair position.
        :return: The tuple (position, is_missing, samples).
        """
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        bpp, fp_data = struct.unpack("QQ", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES))
        flags = (bpp & BpPosFlags.MASK.value)
        is_missing = (0 != (flags & BpPosFlags.IS_MISSING.value))
        is_sparse = (0 != (flags & BpPosFlags.SPARSE.value))
        position = bpp & ~BpPosFlags.MASK.value
        self.file_obj.seek(fp_data)
        if is_sparse:
            sample_list = _read_u32_list(self.file_obj)
        else:
            sample_list = []
            byte_count = _div_round_up(self.num_samples, 8)
            data = self.file_obj.read(byte_count)
            for i in range(byte_count):
                samples_for_bv(data, i, sample_list)
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
        assert BitVector is not None, "Could not import BitVector; try 'pip install BitVector'"
        self.file_obj.seek(self._get_var_idx_offset(variant_idx))
        bpp, fp_data = struct.unpack("QQ", self.file_obj.read(IGDConstants.INDEX_ENTRY_BYTES))
        flags = (bpp & BpPosFlags.MASK.value)
        is_missing = (0 != (flags & BpPosFlags.IS_MISSING.value))
        is_sparse = (0 != (flags & BpPosFlags.SPARSE.value))
        position = bpp & ~BpPosFlags.MASK.value
        self.file_obj.seek(fp_data)
        byte_count = _div_round_up(self.num_samples, 8)
        if is_sparse:
            as_list = _read_u32_list(self.file_obj)
            sample_bv = BitVector.BitVector(size=byte_count*8)
            for sample_idx in as_list:
                sample_bv[sample_idx] = 1
        else:
            data = self.file_obj.read(byte_count)
            sample_bv = BitVector.BitVector(rawbytes=data)
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
            assert count == self._num_idv, f"Malformed file: {count} individual labels but expected {self._num_idv}"
            for _ in range(count):
                result.append(_read_string(self._version, self.file_obj))
        return result

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
            assert count == self._num_var, f"Malformed file: {count} variant labels but expected {self._num_var}"
            for _ in range(count):
                result.append(_read_string(self._version, self.file_obj))
        return result


@dataclass
class IGDHeader:
    magic: int
    version: int
    ploidy: int
    sparse_threshold: int
    num_variants: int
    num_individuals: int
    flags: int
    fp_index: int
    fp_variants: int
    fp_individualids: int
    fp_variantids: int

    def pack(self) -> bytes:
        return struct.pack(IGDConstants.HEADER_FORMAT,
            self.magic, self.version, self.ploidy, self.sparse_threshold, self.num_variants,
            self.num_individuals, self.flags, self.fp_index, self.fp_variants,
            self.fp_individualids, self.fp_variantids, 0, 0, 0, 0, 0, 0)


def _write_u64(file_obj, value):
    file_obj.write(struct.pack("Q", value))

def _write_u32(file_obj, value):
    file_obj.write(struct.pack("I", value))

def _write_str(file_obj, string):
    _write_u32(file_obj, len(string))
    file_obj.write(string.encode("utf-8"))


class IGDWriter:
    def __init__(self,
                 out_stream: BinaryIO,
                 individuals: int,
                 ploidy: int = 2,
                 phased: bool = True,
                 source: str = "",
                 description: str = "",
                 sparse_threshold: int = 32):
        self.out = out_stream
        self.header = IGDHeader(
            magic=IGDConstants.HEADER_MAGIC, version=IGDConstants.SUPPORTED_FILE_VERSION,
            ploidy=ploidy, sparse_threshold=sparse_threshold, num_variants=0,
            num_individuals=individuals,
            flags=0 if not phased else IGDConstants.FLAG_IS_PHASED,
            fp_index=0, fp_variants=0, fp_individualids=0, fp_variantids=0)
        self.source = source
        self.description = description
        self.ref_alleles = []
        self.alt_alleles = []
        self.index = []
        self.num_samples = individuals * ploidy
        self.should_be_sparse = self.num_samples / sparse_threshold
        # The position of the last variant that we wrote. These have to be in order.
        self.last_var_position = 0

    def write_header(self):
        assert self.out.tell() == 0, "Writing header to wrong location"
        self.out.write(self.header.pack())
        _write_str(self.out, self.source)
        _write_str(self.out, self.description)
    
    @staticmethod
    def _make_index_entry(position: int,
                          is_missing: bool,
                          is_sparse: bool,
                          filepos: int) -> bytes:
        flags = 0
        if is_sparse:
            flags |= BpPosFlags.SPARSE.value
        if is_missing:
            flags |= BpPosFlags.IS_MISSING.value
        encoded_pos = position | flags
        return struct.pack("QQ", encoded_pos, filepos)

    def write_variant(self,
                      position: int,
                      ref_allele: str,
                      alt_allele: str,
                      samples: List[int],
                      is_missing: bool = False):
        assert position >= self.last_var_position, "Out of order variant (must be written in ascending order)"
        self.last_var_position = position
        self.ref_alleles.append(ref_allele)
        self.alt_alleles.append(alt_allele)
        is_sparse = len(samples) <= self.should_be_sparse
        filepos = self.out.tell()
        self.index.append(
            self._make_index_entry(position, is_missing, is_sparse, filepos))
        if is_sparse:
            _write_u32(self.out, len(samples))
            prev_sample = -1
            for sample_idx in samples:
                assert prev_sample < sample_idx, "Sample indexes must be in ascending order"
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
            self.out.write(struct.pack("B"*num_bytes, *data))
        self.header.num_variants += 1
    
    def write_index(self):
        assert len(self.index) == self.header.num_variants
        self.header.fp_index = self.out.tell()
        for item in self.index:
            self.out.write(item)

    def write_variant_info(self):
        assert len(self.ref_alleles) == self.header.num_variants
        assert len(self.alt_alleles) == self.header.num_variants
        self.header.fp_variants = self.out.tell()
        for i in range(self.header.num_variants):
            _write_str(self.out, self.ref_alleles[i])
            _write_str(self.out, self.alt_alleles[i])

    def write_individual_ids(self, labels: List[str]):
        if len(labels) == 0:
            self.header.fp_individualids = 0
        else:
            assert len(labels) == self.header.num_individuals
            self.header.fp_individualids = self.out.tell()
            _write_u64(self.out, len(labels))
            for label in labels:
                _write_str(self.out, label)

    def write_variant_ids(self, labels: List[str]):
        if len(labels) == 0:
            self.header.fp_variantids = 0
        else:
            assert len(labels) == self.header.num_variants
            self.header.fp_variantids = self.out.tell()
            _write_u64(self.out, len(labels))
            for label in labels:
                _write_str(self.out, label)


class IGDTransformer:
    def __init__(self,
                 in_stream: BinaryIO,
                 out_stream: BinaryIO,
                 use_bitvectors: bool = False):
        self.reader = IGDReader(in_stream)
        self.writer = IGDWriter(out_stream, self.reader.num_individuals, self.reader.ploidy,
                                self.reader.is_phased, self.reader.source, self.reader.description)
        self.use_bitvectors = use_bitvectors

    def transform(self):
        self.writer.write_header()
        for i in range(self.reader.num_variants):
            if self.use_bitvectors:
                position, is_missing, samples = self.reader.get_samples_bv(i)
            else:
                position, is_missing, samples = self.reader.get_samples(i)
            samples = self.modify_samples(position, is_missing, samples)
            if samples is not None:
                # If the user is using bitvectors convert back to a sample list before writing.
                if self.use_bitvectors:
                    samples_list = []
                    for s in range(len(samples)):
                        if samples[s]:
                            samples_list.append(s)
                    samples = samples_list
                self.writer.write_variant(position,
                                          self.reader.get_ref_allele(i), 
                                          self.reader.get_alt_allele(i), 
                                          samples)
        self.writer.write_index()
        self.writer.write_variant_info()
        self.writer.write_individual_ids(self.reader.get_individual_ids())
        self.writer.write_variant_ids(self.reader.get_variant_ids())
        self.writer.out.seek(0)
        self.writer.write_header()

    def modify_samples(self, position: int, is_missing: bool, samples: Union[BVType, List[int]]) -> bool:
        raise NotImplementedError("Derived class must implement modify_samples()")