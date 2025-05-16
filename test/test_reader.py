from pyigd import IGDReader, IGDConstants, BpPosFlags
import unittest
import struct
import tempfile
import os
import random

TEST_SOURCE = "SOURCE"
TEST_DESCRIPTION = "DESCRIPTION"


def make_header(
    magic=IGDConstants.HEADER_MAGIC,
    version=IGDConstants.SUPPORTED_FILE_VERSION,
    ploidy=2,
    variants=0,
    individuals=0,
    fp_idx=0,
    fp_vars=0,
    fp_indv=0,
    phased=True,
):
    return struct.pack(
        IGDConstants.HEADER_FORMAT,
        magic,
        version,
        ploidy,
        0,
        variants,
        individuals,
        0,
        IGDConstants.FLAG_IS_PHASED if phased else 0,
        fp_idx,
        fp_vars,
        fp_indv,
        0,
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


def _write_str(ver3, file_obj, string):
    if ver3:
        _write_u64(file_obj, len(string))
    else:
        _write_u32(file_obj, len(string))
    file_obj.write(string.encode("utf-8"))


class IGDTestFile(tempfile.TemporaryDirectory):
    FILENAME = "test.igd"

    def __init__(self, header, variants=[], indiv_ids=[], ver3=False):
        super().__init__()
        self.header = header
        self.variants = variants
        self.ver3 = ver3
        self.indiv_ids = indiv_ids

    def __enter__(self, *args):
        tmpdir = super().__enter__(*args)

        outfile = os.path.join(tmpdir, self.FILENAME)
        with open(outfile, "wb") as f:
            # Write the fixed-length header
            f.write(self.header)

            # Write the two strings immediately following header (required)
            _write_str(self.ver3, f, TEST_SOURCE)
            _write_str(self.ver3, f, TEST_DESCRIPTION)

            # Write the genotype data first.
            file_positions = []
            for _, _, sample_list, _ in self.variants:
                file_positions.append(f.tell())
                _write_u32(f, len(sample_list))
                for sample_idx in sample_list:
                    _write_u32(f, sample_idx)

            # Then write the index.
            fp_idx = f.tell()
            for (position, is_missing, _, _), fp in zip(self.variants, file_positions):
                flags = BpPosFlags.SPARSE.value
                if is_missing:
                    flags |= BpPosFlags.IS_MISSING.value
                encoded_pos = position | flags
                _write_u64(f, encoded_pos)
                _write_u64(f, fp)

            # If provided write the individual ids
            if self.indiv_ids:
                fp_indivs = f.tell()
                _write_u64(f, len(self.indiv_ids))
                for ident in self.indiv_ids:
                    _write_str(self.ver3, f, ident)
            else:
                fp_indivs = 0

            # If provided write the variant ids
            var_ids = [ident for _, _, _, ident in self.variants if ident is not None]
            if var_ids:
                fp_varids = f.tell()
                _write_u64(f, len(var_ids))
                for ident in var_ids:
                    _write_str(self.ver3, f, ident)
            else:
                fp_varids = 0

            # Just use A/G as every reference/alternate.
            fp_vars = f.tell()
            for _ in range(len(self.variants)):
                _write_str(self.ver3, f, "A")  # REF
                _write_str(self.ver3, f, "G")  # ALT

            # Gross. If this wasn't test code I wouldn't do this...
            # We're just updating the last few fields of the header to set the file locations.
            f.seek(0 + (6 * 8))
            _write_u64(f, fp_idx)
            _write_u64(f, fp_vars)
            _write_u64(f, fp_indivs)
            _write_u64(f, fp_varids)
            f.flush()

        return tmpdir


class ReaderTests(unittest.TestCase):
    def test_good_header_no_data(self):
        with IGDTestFile(make_header()) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                self.assertEqual(igd_file.description, TEST_DESCRIPTION)
                self.assertEqual(igd_file.source, TEST_SOURCE)
                self.assertEqual(igd_file.num_individuals, 0)
                self.assertEqual(igd_file.num_variants, 0)
                self.assertEqual(igd_file.ploidy, 2)

    def test_good_header_some_data(self):
        I = 10
        V = 50

        def rand_samples():
            num_samples = I * 2
            count = random.randint(0, num_samples - 1)
            result = list(range(num_samples))
            random.shuffle(result)
            return sorted(result[:count])

        variants = list(
            sorted(
                [
                    (random.randint(0, 100_000), False, rand_samples(), None)
                    for i in range(V)
                ],
                key=lambda v: v[0],
            )
        )

        with IGDTestFile(make_header(individuals=I, variants=V), variants) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                self.assertEqual(igd_file.description, TEST_DESCRIPTION)
                self.assertEqual(igd_file.source, TEST_SOURCE)
                self.assertEqual(igd_file.num_individuals, I)
                self.assertEqual(igd_file.num_variants, V)
                self.assertEqual(igd_file.ploidy, 2)

                for i in range(V):
                    pos, m, samples = igd_file.get_samples(i)
                    orig_pos, orig_m, orig_samples, _ = variants[i]
                    self.assertEqual(pos, orig_pos)
                    self.assertEqual(m, orig_m)
                    self.assertEqual(samples, orig_samples)

                self.assertEqual(len(igd_file.get_individual_ids()), 0)

    def test_good_indiv_ids(self):
        I = 10
        with IGDTestFile(
            make_header(individuals=I), indiv_ids=list(map(str, range(I)))
        ) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                self.assertEqual(igd_file.num_individuals, I)
                self.assertEqual(
                    igd_file.get_individual_ids(), list(map(str, range(I)))
                )

    def test_bad_indiv_ids(self):
        I = 10
        with IGDTestFile(
            make_header(individuals=I), indiv_ids=list(map(str, range(I + 1)))
        ) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                with self.assertRaises(AssertionError):
                    igd_file.get_individual_ids()

    def test_bad_var_ids(self):
        V = 25
        variants = [(i, False, [], (None if i > 5 else str(i))) for i in range(V)]
        with IGDTestFile(make_header(variants=V), variants=variants) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                with self.assertRaises(AssertionError):
                    igd_file.get_variant_ids()

    def test_good_var_ids(self):
        V = 21
        variants = [(i, False, [], str(i)) for i in range(V)]
        with IGDTestFile(make_header(variants=V), variants=variants) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                self.assertEqual(igd_file.get_variant_ids(), [v[3] for v in variants])

    def test_lower_bound(self):
        V = 21
        variants = [(i * 10, False, [], str(i)) for i in range(V)]
        with IGDTestFile(make_header(variants=V), variants=variants) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)

                self.assertEqual(0, igd_file.lower_bound_position(0))
                self.assertEqual(
                    igd_file.num_variants, igd_file.lower_bound_position(100_000)
                )
                self.assertEqual(1, igd_file.lower_bound_position(10))
                self.assertEqual(9, igd_file.lower_bound_position(90))
                self.assertEqual(10, igd_file.lower_bound_position(91))
                self.assertEqual(10, igd_file.lower_bound_position(92))
                self.assertEqual(10, igd_file.lower_bound_position(93))
                self.assertEqual(10, igd_file.lower_bound_position(99))
                self.assertEqual(10, igd_file.lower_bound_position(100))
                self.assertEqual(11, igd_file.lower_bound_position(101))


class CompatTests(unittest.TestCase):
    def test_bad_header_ver(self):
        with IGDTestFile(make_header(version=1)) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                with self.assertRaises(AssertionError):
                    _ = IGDReader(f)

    def test_bad_header_magic(self):
        with IGDTestFile(make_header(magic=0xBAADF00D)) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                with self.assertRaises(AssertionError):
                    _ = IGDReader(f)

    def test_v3_compatible(self):
        I = 10
        V = 50
        with IGDTestFile(
            make_header(version=3, individuals=I, variants=V), ver3=True
        ) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)

                self.assertEqual(igd_file.description, TEST_DESCRIPTION)
                self.assertEqual(igd_file.source, TEST_SOURCE)
                self.assertEqual(igd_file.num_individuals, I)
                self.assertEqual(igd_file.num_variants, V)
                self.assertEqual(igd_file.ploidy, 2)


if __name__ == "__main__":
    unittest.main()
