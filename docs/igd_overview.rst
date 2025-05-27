Overview of IGD
===============

An IGD file contains variant positions, identifiers (e.g., "ID" field from VCF), genotype data for a list of samples,
and identifiers for the individuals corresponding to those samples. The genotype data is referenced by *sample index*.
The sample index ``i`` ranges from ``0`` to ``N-1``, where ``N`` is the number of haplotypes. Sample indices are grouped
by individual, so given a ploidy ``P``, every ``P`` consecutive sample indices will be for the same individual.

IGD files are read with :py:class:`pyigd.IGDReader`. They can be written with :py:class:`pyigd.IGDWriter`. There is a helpful
class :py:class:`pyigd.IGDTransformer` that lets you modify an existing IGD file (creating a copy of it) in very few lines
of code. See the `examples <examples.html>`_ or the `API documentation <pyigd.html>`_ for details.

Variants
--------

Variants in IGD are stored bi-allelic, with each IGD variant storing the sample list for a single alternate allele. So if a
site is multi-allelic with ``K`` alternate alleles, there will be ``K`` variants in the IGD file, each with the same position
and reference allele, but different alternate alleles and different sample lists.

The IGD file contains a position table, which can be scanned much faster than scanning all of the genotype data. This table
contains the position and flags for each variant, and points to the location of the genotype data. The index can be searched
linearly (:py:meth:`pyigd.IGDReader.get_position_and_flags`) or via binary search (:py:meth:`pyigd.IGDReader.lower_bound_position`).

Each variant (optionally) has an identifier, which can be looked up via the index of the variant in the array returned by
:py:meth:`pyigd.IGDReader.get_variant_ids`.

Variant Indices
~~~~~~~~~~~~~~~

Everything variant related (genotype data, variant identifiers, reference alleles, alternate alleles) in an IGD file is
looked up by the variant index. Given ``V`` variants in the file, each variant is indexed by a number between ``0`` and ``V-1``.
The variants are ordered according to ascending base-pair position on the chromosome (IGD does not sort them - the file must
be constructed with this order). So the first polymorphic position on the chromosome (for the given dataset) is given
by ``IGDReader.get_position_and_flags(0)`` and the last position is ``IGDReader.get_position_and_flags(IGDReader.num_variants - 1)``.

Individuals
-----------

Given ``N`` haplotype samples and a ploidy of ``P``, there will always be ``N/P`` individuals (use :py:attr:`pyigd.IGDReader.num_individuals`).
The (optional) identifiers for these individuals can be retrieved via :py:meth:`pyigd.IGDReader.get_individual_ids`.

Genotype Data
-------------

The genotype data is retrieved as a list of sample indexes (the ones which contain the alternate allele).

Phased
~~~~~~

For phased data, each sample index corresponds to a haplotype sample, not an individual. The :py:meth:`pyigd.IGDReader.get_samples`
method called on a variant index `i` will return the list of samples that have the alternate allele, which can be retrieved
via :py:meth:`pyigd.IGDReader.get_alt_allele`.

Unphased
~~~~~~~~

For *unphased* data, each sample index corresponds to an individual, not a haplotype. Each variant has an additional piece
of information associated with it: the number of copies of the alternate allele that the individuals have. The number of
copies is between ``1`` and ``P`` (the ploidy), and is obtained by :py:meth:`pyigd.IGDReader.get_position_flags_copies
(returns ``num_copies`` as third item in the returned tuple). For example, with diploid individuals the homozygous
individuals are returned when ``num_copies=2`` and the heterozygous individuals are returns when ``num_copies=1``. When an
individual is homozygous in the reference allele, they will not be in any sample list (homozygous for reference is the implicit
case).  The :py:meth:`pyigd.IGDReader.get_alt_allele` function is still used to retrieve the corresponding sample lists.

Metadata
--------

The IGD file itself does not contain metadata (beyond variant and individual identifiers). However, ``igdtools`` supports
exporting variant-based metadata to files that can be loaded with `numpy.loadtxt <https://numpy.org/doc/2.2/reference/generated/numpy.loadtxt.html>`_.
Matrix-based metadata (i.e., for VCF this means FORMAT fields other than GT) is not supported: if you need per-variant-per-sample metadata, then there
is probably no reason to use IGD (you need a non-sparse representation like VCF/BCF, since your metadata is non-sparse).

There are two ways to export this metadata:

1. During VCF-to-IGD conversion: ``igdtools in.vcf.gz -o out.igd -e all``

2. Only export metadata from a VCF: ``igdtools in.vcf.gz -e all``

The metadata is stored as a file per metadata item type. The supported fields are CHROM, QUAL, FILTER, and INFO. For INFO, each
key gets its own file.  All metadata files are a single entry (line) per variant in the resulting IGD file (i.e., "expanded" variants).
The files are stored in a ``<prefix>.meta/`` directory, where the prefix is determined by the output file (for conversion) or input
file (when you are just exporting metadata).

The first line of a metadata file is a comment that has information about the metadata. When loaded with ``numpy.loadtxt()``, the size of
the array is exactly :py:meth:`pyigd.IGDReader.num_variants` in length, and if you index variant ``i`` in the IGD file you can get its metadata by
looking at element ``i`` of the corresponding metadata array.

When a metadata value is not provided for a particular variant, a default value is used based on the Type field in the VCF metadata:

* Integer: ``0``
* Float: ``NaN``
* String: ``.``

Below is some example Python code for loading metadata files, here we assume we had that INFO field ``AC`` (allele counts) in the original
VCF file that we converted and exported metadata for. 

::

    import pyigd
    import numpy

    with open("test.igd", "rb") as figd, open("test.meta.info.AC.txt") as fmeta:
        igd_file = pyigd.IGDReader(figd)
        ac_meta_data = numpy.loadtxt(fmeta, dtype=int) 
        for i in range(igd_file.num_variants):
            position, flags = igd_file.get_position_and_flags(i)
            ac_value = ac_meta_data[i]
            print(f"Variant={i}, Position={position}, AC={ac_value}")


When an IGD file has been modified/filtered after the metadata was exported, you'll need to match up the variant identifiers between the
two files. When metadata is exported, there is always a ``variants.txt`` file that is the variant identifiers associated with the metadata
at the time of export, so you just need to match those identifiers up with whatever identifiers your IGD file contains, to find the correct
metadata rows. See the `examples page <examples.html>`_ for an example that uses `numpy.intersect1d <https://numpy.org/doc/stable/reference/generated/numpy.intersect1d.html>`_
to do this efficiently.