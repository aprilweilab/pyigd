"""
Module for reading Indexable Genotype Data (IGD) files.
"""
from contextlib import AbstractContextManager
import struct
from enum import Enum
from typing import BinaryIO, List, Tuple
# Optional import. BitVector is only needed when using the APIs that return them.
try:
    import BitVector  # type: ignore
except ImportError:
    BitVector = None


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


def _read_string(file_obj: BinaryIO) -> str:
    length = _read_uint64(file_obj)
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


class IGDFile(AbstractContextManager):
    NUM_HEADER_BYTES = 128
    INDEX_ENTRY_BYTES = 16

    HEADER_FORMAT = "QQIIQQQQQQQQQQQQQ"
    HEADER_MAGIC = 0x3a0c6fd7945a3481
    SUPPORTED_FILE_VERSION = 3

    def __init__(self, filename: str):
        """
        Construct an IGDFile object for reading data from a file.

        :param filename: The path to the file to be opened.
        """
        self.file_obj = open(filename, "rb")

        header_bytes = self.file_obj.read(self.NUM_HEADER_BYTES)
        (magic, self._version, self._ploidy, _, self._num_var, self._num_idv,
         self._flags, self._fp_idx, self._fp_vars, self._fp_ind_ids, _, _, _, _, _,
         _, _) = struct.unpack(self.HEADER_FORMAT, header_bytes)
        assert magic == self.HEADER_MAGIC, "Invalid magic number; not an IGD file"
        assert self._version == self.SUPPORTED_FILE_VERSION, f"Unsupported IGD file format verison {self._version}"

        self._source = _read_string(self.file_obj)
        self._description = _read_string(self.file_obj)

        self._before_first_var = self.file_obj.tell()
        self._all_refs = None
        self._all_alts = None

    def __del__(self):
        if self.file_obj is not None:
            self.file_obj.close()
        self.file_obj = None

    def __exit__(self, *args):
        if self.file_obj is not None:
            self.file_obj.close()
        self.file_obj = None

    def _read_allele_info(self):
        self.file_obj.seek(self._fp_vars)
        self._all_refs = []
        self._all_alts = []
        for i in range(self.num_variants):
            ref = _read_string(self.file_obj)
            self._all_refs.append(ref)
            alt = _read_string(self.file_obj)
            self._all_alts.append(alt)

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
        return self._fp_idx + (variant_idx * self.INDEX_ENTRY_BYTES)

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
        bpp, _ = struct.unpack("QQ", self.file_obj.read(self.INDEX_ENTRY_BYTES))
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
        bpp, fp_data = struct.unpack("QQ", self.file_obj.read(self.INDEX_ENTRY_BYTES))
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
        bpp, fp_data = struct.unpack("QQ", self.file_obj.read(self.INDEX_ENTRY_BYTES))
        flags = (bpp & BpPosFlags.MASK.value)
        is_missing = (0 != (flags & BpPosFlags.IS_MISSING.value))
        is_sparse = (0 != (flags & BpPosFlags.SPARSE.value))
        position = bpp & ~BpPosFlags.MASK.value
        self.file_obj.seek(fp_data)
        if is_sparse:
            as_list = _read_u32_list(self.file_obj)
            max_index = as_list[-1] if as_list else 0
            sample_bv = BitVector.BitVector(size=max_index+1)
            for sample_idx in as_list:
                sample_bv[sample_idx] = 1
        else:
            byte_count = _div_round_up(self.num_samples, 8)
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

    def get_individual_ids(self):
        raise NotImplementedError("Retrieving individual IDs from IGD file not yet supported")
