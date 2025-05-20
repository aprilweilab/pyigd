from pyigd import IGDReader, IGDWriter
from pyigd.extra import collect_next_site, get_inverse_sample_list
import unittest
import tempfile
import os

TEST_SOURCE = "SOURCE"
TEST_DESCRIPTION = "DESCRIPTION"


# 20 individuals
class IGDTestFile20(tempfile.TemporaryDirectory):
    FILENAME = "test.igd"
    INDIVS = 20

    def __init__(self, positions, refs, alts, samples, num_copies=[]):
        super().__init__()
        assert len(positions) == len(refs)
        assert len(alts) == len(refs)
        assert len(samples) == len(alts)
        if not num_copies:
            num_copies = [0] * len(alts)
        assert len(num_copies) == len(alts)
        self.variants = list(zip(positions, refs, alts, samples, num_copies))

    def __enter__(self, *args):
        tmpdir = super().__enter__(*args)

        outfile = os.path.join(tmpdir, self.FILENAME)
        unphased = any(map(lambda v: v[4] > 0, self.variants))
        with open(outfile, "wb") as f:
            writer = IGDWriter(f, individuals=self.INDIVS, phased=not unphased)
            writer.write_header()
            for pos, ref, alt, samples, num_copies in self.variants:
                is_missing = False
                if len(alt) == 0:
                    is_missing = True
                writer.write_variant(pos, ref, alt, samples, is_missing, num_copies)
            writer.write_index()
            writer.write_variant_info()
            writer.out.seek(0)
            writer.write_header()

        return tmpdir


class ReaderTests(unittest.TestCase):
    def test_site_iterate_phased(self):
        positions = [1, 1, 2, 2, 3, 3, 3]
        refs = ["A", "A", "T", "T", "C", "C", "C"]
        alts = ["G", "C", "G", "A", "T", "A", "G"]
        samples = [
            [1, 5],
            [2, 6],
            list(range(10)),
            list(range(10, 20)),
            list(range(0, 5)),
            list(range(10, 15)),
            list(range(25, 30)),
        ]
        with IGDTestFile20(positions, refs, alts, samples) as tmpdir:
            result = []
            filename = os.path.join(tmpdir, IGDTestFile20.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                i = 0
                while i < igd_file.num_variants:
                    site_variants = collect_next_site(igd_file, i)
                    i = site_variants[-1] + 1
                    result.append(site_variants)
                self.assertEqual(
                    result,
                    [
                        [0, 1],
                        [2, 3],
                        [4, 5, 6],
                    ],
                )
                all_ref_samples = []
                for indices in result:
                    ref_samples = get_inverse_sample_list(igd_file, indices)
                    assert len(ref_samples) == 1  # Phased data
                    all_ref_samples.append(ref_samples[0])

                self.assertEqual(
                    all_ref_samples,
                    [
                        [0, 3, 4] + list(range(7, 40)),
                        list(range(20, 40)),
                        list(range(5, 10)) + list(range(15, 25)) + list(range(30, 40)),
                    ],
                )

    def test_site_iterate_unphased(self):
        positions = [1, 1, 1, 2, 2, 3, 3, 3]
        refs = ["A", "A", "A", "T", "T", "C", "C", "C"]
        alts = ["G", "", "C", "A", "A", "T", "A", "T"]
        num_copies = [1, 2, 1, 1, 2, 1, 1, 2]
        samples = [
            [1, 5],
            [2, 6],
            [5, 12, 13],
            [1, 2, 3, 4, 5, 6, 7],
            [18],
            [18, 19],
            [11, 18, 19, 15],
            [1, 2, 3, 4],
        ]
        with IGDTestFile20(positions, refs, alts, samples, num_copies) as tmpdir:
            result = []
            filename = os.path.join(tmpdir, IGDTestFile20.FILENAME)
            with open(filename, "rb") as f:
                igd_file = IGDReader(f)
                i = 0
                while i < igd_file.num_variants:
                    site_variants = collect_next_site(igd_file, i)
                    i = site_variants[-1] + 1
                    result.append(site_variants)
                self.assertEqual(
                    result,
                    [
                        [0, 1, 2],
                        [3, 4],
                        [5, 6, 7],
                    ],
                )

                ref_hets = []
                ref_homs = []
                for indices in result:
                    ref_samples = get_inverse_sample_list(igd_file, indices)
                    # These three should be the same.
                    assert len(ref_samples) == 2
                    assert igd_file.ploidy == len(ref_samples)

                    ref_hets.append(ref_samples[0])
                    ref_homs.append(ref_samples[1])

                self.assertEqual(
                    ref_hets,
                    [
                        [1, 12, 13],  # First site
                        [1, 2, 3, 4, 5, 6, 7],  # Second site
                        [11, 15],  # Third site
                    ],
                )
                self.assertEqual(
                    ref_homs,
                    [
                        [
                            0,
                            3,
                            4,
                            7,
                            8,
                            9,
                            10,
                            11,
                            14,
                            15,
                            16,
                            17,
                            18,
                            19,
                        ],  # First site
                        [0, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19],  # Second site
                        [0, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17],  # Third site
                    ],
                )
