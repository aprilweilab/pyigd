from pyigd import IGDReader, IGDWriter, IGDTransformer
import unittest
import tempfile
import os

TEST_SOURCE = "SOURCE"
TEST_DESCRIPTION = "DESCRIPTION"

class CanonicalTestIGD:
    def __init__(self):
        self.variants = [
            (1, "A", "G", list(range(500, 1500))),
            (1, "C", "G", list(range(0, 2000))),
            (1, "T", "REALLYLONG", [1,2]),
            (1, "A", "T", []),
            (1, "G", "G", [0, 9, 19, 30, 42, 55]),
        ]
    
    def write(self, filename: str):
        with open(filename, "wb") as f:
            w = IGDWriter(f, 2000, ploidy=1, phased=False, source="TEST", description="TESTD")
            w.write_header()
            for pos, ref, alt, samples in self.variants:
                w.write_variant(pos, ref, alt, samples)
            w.write_index()
            w.write_variant_info()
            w.write_individual_ids([f"ID{i}" for i in range(2000)])
            w.write_variant_ids([f"VAR{i}" for i in range(len(self.variants))])
            f.seek(0)
            w.write_header()
            
    def readAndVerify(self, filename: str, test_case: unittest.TestCase):
        with open(filename, "rb") as f:
            reader = IGDReader(f)
            test_case.assertEqual(reader.num_individuals, 2000)
            test_case.assertEqual(reader.num_variants, len(self.variants))
            for i, (pos, ref, alt, samples) in enumerate(self.variants):
                test_case.assertEqual(reader.get_alt_allele(i), alt)
                test_case.assertEqual(reader.get_ref_allele(i), ref)
                rpos, is_missing, rsamples = reader.get_samples(i)
                test_case.assertEqual(rpos, pos)
                test_case.assertFalse(is_missing)
                test_case.assertEqual(rsamples, samples)
            test_case.assertEqual(reader.description, "TESTD")
            test_case.assertEqual(reader.source, "TEST")


class WriterTests(unittest.TestCase):
    def test_write_read_simple(self):
        """
        Write a simple IGD file verify reading.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            fn = os.path.join(tmpdir, "tmp.igd")
            with open(fn, "wb") as f:
                w = IGDWriter(f, 100)
                w.write_header()
                w.write_variant(100, "A", "G", [1, 22, 99, 101])
                w.write_variant(101, "A", "C", list(range(200)))
                w.write_index()
                w.write_variant_info()
                w.write_individual_ids([])
                w.write_variant_ids([])
                f.seek(0)
                w.write_header()
            
            with open(fn, "rb") as f:
                reader = IGDReader(f)
                self.assertEqual(reader.num_individuals, 100)
                self.assertEqual(reader.num_variants, 2)
                self.assertEqual(reader.get_alt_allele(0), "G")
                self.assertEqual(reader.get_alt_allele(1), "C")
                self.assertEqual(reader.get_samples(0), (100, False, [1, 22, 99, 101]))
                self.assertEqual(reader.get_samples(1), (101, False, list(range(200))))

    def test_write_read_loop(self):
        """
        Write an IGD file with as many interesting features as possible and verify that reading
        recovers all of that data.
        """
        test_data = CanonicalTestIGD()
        with tempfile.TemporaryDirectory() as tmpdir:
            fn = os.path.join(tmpdir, "test.igd")
            test_data.write(fn)
            test_data.readAndVerify(fn, self)

    def test_out_of_order(self):
        """
        Writing should fail if variants are not in ascending positional order.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            fn = os.path.join(tmpdir, "tmp.igd")
            with open(fn, "wb") as f:
                w = IGDWriter(f, 200)
                w.write_variant(100, "A", "G", [1])
                with self.assertRaises(AssertionError) as context:
                    w.write_variant(99, "A", "G", [1])

    def test_bad_samples(self):
        """
        Writing should fail if a bad sample ID is given.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            fn = os.path.join(tmpdir, "tmp.igd")
            with open(fn, "wb") as f:
                w = IGDWriter(f, 200)
                with self.assertRaises(AssertionError) as context:
                    w.write_variant(100, "A", "G", [99, 98])

class TransformerTests(unittest.TestCase):
    def test_copy(self):
        """
        Copy an IGD exactly using IGDTransformer.
        """
        class MyXformer(IGDTransformer):
            def modify_samples(self, position, is_missing, samples):
                return samples

        test_data = CanonicalTestIGD()
        with tempfile.TemporaryDirectory() as tmpdir:
            in_file = os.path.join(tmpdir, "created.igd")
            out_file = os.path.join(tmpdir, "copied.igd")
            out_file_bv = os.path.join(tmpdir, "copied_bv.igd")

            test_data.write(in_file)
            with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
                xformer = MyXformer(fin, fout)
                xformer.transform()
            test_data.readAndVerify(out_file, self)

            with open(in_file, "rb") as fin, open(out_file_bv, "wb") as fout:
                xformer = MyXformer(fin, fout, use_bitvectors=True)
                xformer.transform()
            test_data.readAndVerify(out_file_bv, self)

    def test_simple_mod(self):
        """
        Copy an IGD dropping the first sample of every list, using IGDTransformer.
        """
        class MyXformer(IGDTransformer):
            def modify_samples(self, position, is_missing, samples):
                return samples[1:]

        test_data = CanonicalTestIGD()
        with tempfile.TemporaryDirectory() as tmpdir:
            in_file = os.path.join(tmpdir, "created.igd")
            out_file = os.path.join(tmpdir, "copied.igd")

            test_data.write(in_file)
            with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
                xformer = MyXformer(fin, fout)
                xformer.transform()
            # Before verifying, remove the first sample from our test data as well.
            for i in range(len(test_data.variants)):
                test_sample_list = test_data.variants[i][3]
                if test_sample_list:
                    del test_sample_list[0]
            # The modified file should match our test data now
            test_data.readAndVerify(out_file, self)
            # The original file should _not_ match out test data now
            with self.assertRaises(AssertionError) as context:
                test_data.readAndVerify(in_file, self)

    def test_simple_mod_bv(self):
        """
        Copy an IGD dropping the first sample of every list, using IGDTransformer.
        """
        class MyXformer(IGDTransformer):
            def modify_samples(self, position, is_missing, samples):
                for i in range(len(samples)):
                    if samples[i]:
                        samples[i] = 0
                        break
                return samples

        test_data = CanonicalTestIGD()
        with tempfile.TemporaryDirectory() as tmpdir:
            in_file = os.path.join(tmpdir, "created.igd")
            out_file = os.path.join(tmpdir, "copied.igd")

            test_data.write(in_file)
            with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
                xformer = MyXformer(fin, fout, use_bitvectors=True)
                xformer.transform()
            # Before verifying, remove the first sample from our test data as well.
            for i in range(len(test_data.variants)):
                test_sample_list = test_data.variants[i][3]
                if test_sample_list:
                    del test_sample_list[0]
            # The modified file should match our test data now
            test_data.readAndVerify(out_file, self)
            # The original file should _not_ match out test data now
            with self.assertRaises(AssertionError) as context:
                test_data.readAndVerify(in_file, self)



if __name__ == "__main__":
    unittest.main()