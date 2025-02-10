# Example: skip all missing data.
import pyigd
import sys

if len(sys.argv) < 2:
    print("Pass in IGD filename")
    exit(1)

with open(sys.argv[1], "rb") as f:
    igd_file = pyigd.IGDReader(f)
    print(f"Version: {igd_file.version}")
    print(f"Ploidy: {igd_file.ploidy}")
    for variant_index in range(igd_file.num_variants):
        position, flags = igd_file.get_position_and_flags(variant_index)
        if pyigd.flags_is_missing(flags):
            continue
        _, _, sample_list = igd_file.get_samples(variant_index)
        print(f"{position}: {len(sample_list)} samples")
