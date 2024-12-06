# Example demonstrating runs-of-homozygosity using unphased IGD data.
import pyigd
import sys

# 500kbp
THRESHOLD = 500_000

if len(sys.argv) < 2:
    print("Pass in an IGD filename")
    exit(1)

with open(sys.argv[1], "rb") as f:
    igd_file = pyigd.IGDReader(f)
    assert not igd_file.is_phased, "This example only works for unphased data"
    print(f"Version: {igd_file.version}", file=sys.stderr)
    print(f"Ploidy: {igd_file.ploidy}", file=sys.stderr)
    print(f"Variants: {igd_file.num_variants}", file=sys.stderr)
    print(f"Individuals: {igd_file.num_individuals}", file=sys.stderr)
    print(f"Source: {igd_file.source}", file=sys.stderr)
    print(f"Description: {igd_file.description}", file=sys.stderr)
    assert igd_file.num_samples == igd_file.num_individuals

    last_het_site_per_idv = [0 for _ in range(igd_file.num_individuals)]
    print(f"INDIV\tSTART_BP\tEND_BP")
    for variant_index in range(igd_file.num_variants):
        position, flags, num_copies = igd_file.get_position_flags_copies(variant_index)
        # Homozygous occurs when the number of copies of the alt allele are equal to the
        # organism ploidy. We check all heterozygous cases as they "break" the ROH, so
        # if we break an ROH that exceeds the threshold we emit it.
        if num_copies < igd_file.ploidy:
            _, is_missing, sample_list = igd_file.get_samples(variant_index)
            assert not is_missing, "Missing data not handled in this example"
            for indiv in sample_list:
                hom_span = position - last_het_site_per_idv[indiv]
                if hom_span >= THRESHOLD:
                    print(f"{indiv}\t{last_het_site_per_idv[indiv]+1}\t{position-1}")
                last_het_site_per_idv[indiv] = position
    for indiv in range(igd_file.num_individuals):
        hom_span = position - last_het_site_per_idv[indiv]
        if hom_span >= THRESHOLD:
            print(f"{indiv}\t{last_het_site_per_idv[indiv]+1}\t{position-1}")
