import pyigd
import sys

if len(sys.argv) < 3:
    print("Pass in IGD input/output filenames")
    exit(1)


# Our own IGDTransformer. This just deletes the first sample from every variant, but
# in general you can change the sample list any way you want. Return None to delete
# the variant entirely.
class MyXformer(pyigd.IGDTransformer):
    def modify_samples(self, position, is_missing, samples, num_copies):
        return samples[1:]


with open(sys.argv[1], "rb") as fin, open(sys.argv[2], "wb") as fout:
    xformer = MyXformer(fin, fout, use_bitvectors=False)
    xformer.transform()
print(
    f"Copied {sys.argv[1]} to {sys.argv[2]} and deleted the first sample from every variant."
)
