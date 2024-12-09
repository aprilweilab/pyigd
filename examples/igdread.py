import pyigd
import sys

if len(sys.argv) < 2:
    print("Pass in IGD filename")
    exit(1)

with open(sys.argv[1], "rb") as f:
    igd_file = pyigd.IGDReader(f)
    print(f"Version: {igd_file.version}")
    print(f"Ploidy: {igd_file.ploidy}")
    print(f"Variants: {igd_file.num_variants}")
    print(f"Individuals: {igd_file.num_individuals}")
    print(f"Source: {igd_file.source}")
    print(f"Description: {igd_file.description}")
    for variant_index in range(igd_file.num_variants):
        # Approach 1: Get the samples as a list
        print(
            f"REF: {igd_file.get_ref_allele(variant_index)}, ALT: {igd_file.get_alt_allele(variant_index)}"
        )
        position, is_missing, sample_list = igd_file.get_samples(variant_index)
        print((position, is_missing, len(sample_list)))

        # Approach 2: Get the samples as a BitVector object
        # See https://engineering.purdue.edu/kak/dist/BitVector-3.5.0.html
        position, is_missing, bitvect = igd_file.get_samples_bv(variant_index)
        print((position, is_missing, bitvect.count_bits()))
