![Python build and test](https://github.com/aprilweilab/pyigd/actions/workflows/python-package.yml/badge.svg)

# pyigd

PyIGD is a Python-only parser the [Indexable Genotype Data (IGD) format](https://github.com/aprilweilab/picovcf/blob/main/IGD.FORMAT.md).

For a C++ library that supports creating and parsing IGD, see [picovcf](https://github.com/aprilweilab/picovcf) (which also supports VCF -> IGD conversion).

## Installation

Clone the code and then either install for development:
```
pip install -e pyigd/
```

or build and install via the wheel:
```
cd pyigd/ && python setup.py bdist_wheel
pip install --force-reinstall dist/*.whl
```

## Usage

The `pyigd.IGDFile` class is a context manager, so it is recommended that you use it via the `with` statement. Below is an example script that loads an IGD file, prints out some meta-data, and then iterates the genotype data for all variants.

```
import pyigd
import sys

if len(sys.argv) < 2:
    print("Pass in IGD filename")
    exit(1)

with pyigd.IGDFile(sys.argv[1]) as igd_file:
    print(f"Version: {igd_file.version}")
    print(f"Ploidy: {igd_file.ploidy}")
    print(f"Variants: {igd_file.num_variants}")
    print(f"Individuals: {igd_file.num_individuals}")
    print(f"Source: {igd_file.source}")
    print(f"Description: {igd_file.description}")
    for variant_index in range(igd_file.num_variants):
        # Approach 1: Get the samples as a list
        print(f"REF: {igd_file.get_ref_allele(variant_index)}, ALT: {igd_file.get_alt_allele(variant_index)}")
        position, is_missing, sample_list = igd_file.get_samples(variant_index)
        print( (position, is_missing, len(sample_list))  )

        # Approach 2: Get the samples as a BitVector object
        # See https://engineering.purdue.edu/kak/dist/BitVector-3.5.0.html
        position, is_missing, bitvect = igd_file.get_samples_bv(variant_index)
        print( (position, is_missing, bitvect.count_bits())  )
```

IGD can be highly performant for a few reasons:
1. It stores sparse data sparsely. Low-frequency variants are stored as sample lists. Medium/high frequency variants are stored as bit vectors.
2. It is indexable (you can jump directly to data for the `ith` variant). Since the index is stored in its own section of the file, scanning the index is extremely fast. So only looking at variants for a particular range of the genome is very fast (in this case you would use `pyigd.IGDFile.get_position_and_flags()` to find the first variant index within the range, and then use `pyigd.IGDFile.get_samples()` after that).
3. The genotype data is stored in one of two very simple binary formats. This makes parsing fast, and the compact nature of the file makes reading from disk/memory fast as well.

## How do I use IGD in my project?

* Clone [picovcf](https://github.com/aprilweilab/picovcf) and follow the instructions in its [README](https://github.com/aprilweilab/picovcf/blob/main/README.md) to build the example tools for that library.
  * If you want to be able to convert `.vcf.gz` (compressed VCF) to IGD, make sure you build with `-DENABLE_VCF_GZ=ON`
* One of the built tools will be `vcfconf`, which converts from VCF to IGD. Run `vcfconv <vcf file> <igd file>` to convert your data to IGD.
* Do one of the following:
  * If your project is C++, copy [picovcf.hpp](https://github.com/aprilweilab/picovcf/blob/main/picovcf.hpp) into your project, `#include` it somewhere and then use according to the [documentation](https://picovcf.readthedocs.io/en/latest/)
  * If your project is Python, clone [pyigd](https://github.com/aprilweilab/pyigd/) and install it per the [README instructions](https://github.com/aprilweilab/pyigd/blob/main/README.md).
