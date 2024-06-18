from pyigd import IGDFile, BpPosFlags
import unittest
import struct
import tempfile
import os
import random

TEST_SOURCE = "SOURCE"
TEST_DESCRIPTION = "DESCRIPTION"

def make_header(magic=IGDFile.HEADER_MAGIC,
                version=IGDFile.SUPPORTED_FILE_VERSION,
                ploidy=2,
                variants=0,
                individuals=0,
                fp_idx=0,
                fp_vars=0,
                fp_indv=0):
    return struct.pack(IGDFile.HEADER_FORMAT,
        magic, version, ploidy, 0, variants, individuals, 0, fp_idx, fp_vars, fp_indv,
        0, 0, 0, 0, 0, 0, 0)

def _write_u64(file_obj, value):
    file_obj.write(struct.pack("Q", value))

def _write_u32(file_obj, value):
    file_obj.write(struct.pack("I", value))

class IGDTestFile(tempfile.TemporaryDirectory):
    FILENAME = "test.igd"

    def __init__(self, header, variants=[]):
        super().__init__()
        self.header = header
        self.variants = variants

    def __enter__(self, *args):
        tmpdir = super().__enter__(*args)

        outfile = os.path.join(tmpdir, self.FILENAME)
        with open(outfile, "wb") as f:
            # Write the fixed-length header
            f.write(self.header)

            # Write the two strings immediately following header (required)
            _write_u64(f, len(TEST_SOURCE))
            f.write(TEST_SOURCE.encode("utf-8"))
            _write_u64(f, len(TEST_DESCRIPTION))
            f.write(TEST_DESCRIPTION.encode("utf-8"))

            # Write the genotype data first.
            file_positions = []
            for _, _, sample_list in self.variants:
                file_positions.append(f.tell())
                _write_u32(f, len(sample_list))
                for sample_idx in sample_list:
                    _write_u32(f, sample_idx)

            # Then write the index.
            fp_idx = f.tell()
            for (position, is_missing, _), fp in zip(self.variants, file_positions):
                flags = BpPosFlags.SPARSE.value
                if is_missing:
                    flags |= BpPosFlags.IS_MISSING.value
                encoded_pos = position | flags
                _write_u64(f, encoded_pos)
                _write_u64(f, fp)

            # Gross. If this wasn't test code I wouldn't do this...
            # We're just updating one field of the header to set the location of the index.
            f.seek(0 + (6*8))
            _write_u64(f, fp_idx)
            f.flush()

        return tmpdir

class BasicTests(unittest.TestCase):
    def test_good_header_no_data(self):
        with IGDTestFile(make_header()) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            igd_file = IGDFile(filename)
            self.assertEqual(igd_file.description, TEST_DESCRIPTION)
            self.assertEqual(igd_file.source, TEST_SOURCE)
            self.assertEqual(igd_file.num_individuals, 0)
            self.assertEqual(igd_file.num_variants, 0)
            self.assertEqual(igd_file.ploidy, 2)

    def test_good_header_some_data(self):
        I = 10
        V = 50

        def rand_samples():
            num_samples = I*2
            count = random.randint(0, num_samples-1)
            result = list(range(num_samples))
            random.shuffle(result)
            return sorted(result[:count])

        variants = list(sorted([(random.randint(0, 100_000), False, rand_samples())
                                for i in range(V)], key=lambda v: v[0]))

        with IGDTestFile(make_header(individuals=I, variants=V), variants) as tmpdir:
            filename = os.path.join(tmpdir, IGDTestFile.FILENAME)
            igd_file = IGDFile(filename)
            self.assertEqual(igd_file.description, TEST_DESCRIPTION)
            self.assertEqual(igd_file.source, TEST_SOURCE)
            self.assertEqual(igd_file.num_individuals, I)
            self.assertEqual(igd_file.num_variants, V)
            self.assertEqual(igd_file.ploidy, 2)

            for i in range(V):
                pos, m, samples = igd_file.get_samples(i)
                orig_pos, orig_m, orig_samples = variants[i]
                self.assertEqual(pos, orig_pos)
                self.assertEqual(m, orig_m)
                self.assertEqual(samples, orig_samples)