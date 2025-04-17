![Python build and test](https://github.com/aprilweilab/pyigd/actions/workflows/python-package.yml/badge.svg)

# pyigd

PyIGD is a Python-only parser the [Indexable Genotype Data (IGD) format](https://github.com/aprilweilab/picovcf/blob/main/IGD.FORMAT.md). We have a short
[preprint paper](https://www.biorxiv.org/content/10.1101/2025.02.05.636549v1.abstract) that describes the format and some of its advantages.

For a C++ library that supports creating and parsing IGD, see [picovcf](https://github.com/aprilweilab/picovcf) (which also supports VCF -> IGD conversion).

## Installation

You can install the latest release of PyIGD from pypi, via `pip install pyigd`.

For development, you can clone the code install it directly from the directory (this will automatically reflect any code changes you make):
```
pip install -e pyigd/
```

or build and install via the wheel:
```
cd pyigd/ && python setup.py bdist_wheel
pip install --force-reinstall dist/*.whl
```

## Usage

The `pyigd.IGDReader` class reads IGD data from a buffer. See [the example script](https://github.com/aprilweilab/pyigd/blob/main/examples/igdread.py) that loads an IGD file, prints out some meta-data, and then iterates the genotype data for all variants. Generally the usage pattern is:
```
with open(filename, "rb") as f:
  igd_reader = pyigd.IGDReader(f)
```

There is also the `pyigd.IGDWriter` class to construct IGD files. Related is `pyigd.IGDTransformer`, which is a way to create a copy of an IGD while modifying its contents. See the IGDTransformer [sample list example](https://github.com/aprilweilab/pyigd/blob/main/examples/xform.py) and [bitvector example](https://github.com/aprilweilab/pyigd/blob/main/examples/xform_bv.py).

IGD can be highly performant for a few reasons:
1. It stores sparse data sparsely. Low-frequency variants are stored as sample lists. Medium/high frequency variants are stored as bit vectors.
2. It is indexable (you can jump directly to data for the `ith` variant). Since the index is stored in its own section of the file, scanning the index is extremely fast. So only looking at variants for a particular range of the genome is very fast (in this case you would use `pyigd.IGDFile.get_position_and_flags()` to find the first variant index within the range, and then use `pyigd.IGDFile.get_samples()` after that).
3. The genotype data is stored in one of two very simple binary formats. This makes parsing fast, and the compact nature of the file makes reading from disk/memory fast as well.

## How do I use IGD in my project?

* Clone [picovcf](https://github.com/aprilweilab/picovcf) and follow the instructions in its [README](https://github.com/aprilweilab/picovcf/blob/main/README.md) to build the tools for that library.
  * If you want to be able to convert `.vcf.gz` (compressed VCF) to IGD, make sure you build with `-DENABLE_VCF_GZ=ON`
* One of the built tools will be `igdtools`, which can converts from VCF to IGD, among other things (such as filtering IGD files).
* Do one of the following:
  * If your project is C++, copy [picovcf.hpp](https://github.com/aprilweilab/picovcf/blob/main/picovcf.hpp) into your project, `#include` it somewhere and then use according to the [documentation](https://picovcf.readthedocs.io/en/latest/)
  * If your project is Python, clone [pyigd](https://github.com/aprilweilab/pyigd/) and install it per the [README instructions](https://github.com/aprilweilab/pyigd/blob/main/README.md).
