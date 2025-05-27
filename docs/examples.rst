Example Usage
=============

IGD Traversal
~~~~~~~~~~~~~

The two core concepts for IGD traversal are the IGD index and the IGD genotype data. The
IGD index is a contiguous file region that only contains the information about position
and flags for each variant, so it can be scanned very quickly. The IGD genotype data has
all the sample information for each variant, and is therefore slower to access. The IGD
Index is accessed via :py:meth:`pyigd.IGDReader.get_position_and_flags` or
:py:meth:`pyigd.IGDReader.get_position_flags_copies`. The genotype data is accessed via
:py:meth:`pyigd.IGDReader.get_samples` or :py:meth:`pyigd.IGDReader.get_samples_bv`.

Traverse while skipping missing data
------------------------------------

Missing data can be identified during traversal of the IGD index or IGD genotype data. We
show the latter first:

::

  import pyigd

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    for variant_index in range(igd_file.num_variants):
      position, is_missing, sample_list = igd_file.get_samples(variant_index)
      if not is_missing:
        pass # Do something with sample_list


If you expect there to be a lot of missing data, the above can be slightly slow, as you are
retrieving the entire sample list for each variant with missing alleles. Here is the more
efficient IGD index version:

::

  import pyigd

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    for variant_index in range(igd_file.num_variants):
      position, flags = igd_file.get_position_and_flags(variant_index)
      if not pyigd.flags_is_missing(flags):
        _, _, sample_list = igd_file.get_samples(variant_index)
        # Do something with sample_list


Variants in a specific range
----------------------------

Similar to the missing data example above, variants in a specific range can be efficiently
found by scanning the IGD index:

::

  import pyigd

  my_range = (5_000_000, 10_000_000)  # 5MBP to 10MBP
  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    for variant_index in range(igd_file.num_variants):
      position, flags = igd_file.get_position_and_flags(variant_index)
      if position >= my_range[0] and position < my_range[1]:
        _, _, sample_list = igd_file.get_samples(variant_index)
        # Do something with sample_list


Or even more efficient is to use :py:meth:`pyigd.IGDReader.lower_bound_position` to binary search for
first position greater-than-or-equal-to the start of the range, and then traverse until the end.

::

  import pyigd

  my_range = (5_000_000, 10_000_000)  # 5MBP to 10MBP
  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    variant_index = igd_file.lower_bound_position(my_range[0])
    while variant_index < igd_file.num_variants:
      position, flags = igd_file.get_position_and_flags(variant_index)
      if position >= my_range[1]:
        break
      _, _, sample_list = igd_file.get_samples(variant_index)
      # Do something with sample_list


Runs of homozygosity
--------------------

Phased and unphased data are both stored in IGD with a row being equal to an alternate allele at a particular
polymorphic site (i.e., a variant). However, unphased data may have multiple rows for each such variant, one
for each number of copies that the unphased individual has. For example, for diploid individuals there may be
up to two rows for each variant: one row for the scenario where an individual has a single copy of the
alternate allele, and one row for the scenario where an individual has two copies of the alternate allele.
This is stored on the IGD index, as the ``num_copies`` field, and can be accessed via the
:py:meth:`pyigd.IGDReader.get_position_flags_copies` method.

Whereas for phased data each sample is a haplotype, for unphased data each sample is an individual. That is,
``igd_file.num_samples == igd_file.num_individuals``.

Here is an example that finds runs of homozygosity beyond some given threshold, by traversing the IGD index:

::

  import pyigd

  THRESHOLD = 500_000 # Only report ROH exceeding 500kbp

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    assert not igd_file.is_phased, "This example is only for unphased data"
    assert igd_file.ploidy == 2, "This example is for diploids only"

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

    # The last ROH may have gone to the end of the chromosome, so we check for those.
    for indiv in range(igd_file.num_individuals):
        hom_span = position - last_het_site_per_idv[indiv]
        if hom_span >= THRESHOLD:
            print(f"{indiv}\t{last_het_site_per_idv[indiv]+1}\t{position-1}")



Print zygosity counts for each variant
--------------------------------------

Here is a similar example, but instead of finding runs of homozygosity we are simply printing out the
zygosity information for each variant.

::

  import pyigd

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    assert not igd_file.is_phased, "This example is only for unphased data"
    assert igd_file.ploidy == 2, "This example is for diploids only"

    print(f"POSITION\tREF\tALT\tAA\tAa\taa")
    for variant_index in range(igd_file.num_variants):
        he_count = 0
        hz_count = 0
        position, is_missing, bitvect = igd_file.get_samples_bv(variant_index)
        for indiv in range(igd_file.num_individuals):
            sample0 = indiv * 2
            sample1 = sample0 + 1
            if bitvect[sample0] and bitvect[sample1]:
                hz_count += 1
            elif bitvect[sample0] or bitvect[sample1]:
                he_count += 1
        ref = igd_file.get_ref_allele(variant_index)
        alt = igd_file.get_alt_allele(variant_index)
        print(f"{position}\t{ref}\t{alt}\t{igd_file.num_individuals - (hz_count+he_count)}\t{he_count}\t{hz_count}")


Find common positions between two IGD files
-------------------------------------------

Here is an example where we compare two IGD files, to see if the polymorphic sites they contain
are the same:

::

  import pyigd

  unique_pos1 = set()
  with open("file1.igd", "rb") as f:
  igd_file = pyigd.IGDReader(f)
  for variant_index in range(igd_file.num_variants):
     position, flags, num_copies = igd_file.get_position_flags_copies(variant_index)
     unique_pos1.add(position)

  unique_pos2 = set()
  with open("file2.igd", "rb") as f:
  igd_file = pyigd.IGDReader(f)
  for variant_index in range(igd_file.num_variants):
     position, flags, num_copies = igd_file.get_position_flags_copies(variant_index)
     unique_pos2.add(position)

  extra1 = unique_pos1 - unique_pos2
  extra2 = unique_pos2 - unique_pos1
  if len(extra1) != 0:
    print(f"File1 has {len(extra1)} extra positions")
  if len(extra2) != 0:
    print(f"File2 has {len(extra2)} extra positions")
  elif len(extra1) == len(extra2):
    print("Positions are identical")


Traverse by site instead of variant
-----------------------------------

Use the :py:meth:`pyigd.extra.collect_next_site` method to more easily traverse the data by site instead of
by variant. This method collects all the variant indices for the next site into a list.

::

  import pyigd
  from pyigd.extra import collect_next_site

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    next_index = 0
    while next_index < igd_file.num_variants:
      variant_indices = collect_next_site(igd_file, next_index)
      # The indices are ordered, so the next time we iterate we start at the last index + 1
      next_index = variant_indices[-1] + 1

      # Count the number of samples that have _any alternate allele_ at the site:
      alt_count = 0
      for index in variant_indices:
        position, is_missing, samples = igd_file.get_samples(index)
        alt_count += len(samples)
      print(f"At site {position} there are {len(variant_indices)} variants and {alt_count} total alternate alleles")


Polarize during traversal
-----------------------------------

This example uses the `pyfaidx <https://pypi.org/project/pyfaidx/>`_ package for reading FASTA files.
We assume you have an input ancestral sequence in FASTA format (like the ones in
`ENSEMBL release v112 <https://ftp.ensembl.org/pub/release-112/fasta/ancestral_alleles/>`_)

You can use :py:meth:`pyigd.extra.collect_next_site` method to traverse the data by site and then
:py:meth:`pyigd.extra.get_inverse_sample_list` to flip the reference allele (if needed).

::

  import pyigd
  import pyfaidx
  from pyigd.extra import collect_next_site, get_inverse_sample_list

  fasta_reader = pyfaidx.Fasta("ancestral.fa")
  assert len(fasta_reader.values()) == 1
  # Assumes FASTA is 1-based, but we need 0-based indexing.
  ancestral_str = ("X" + str(list(fasta_reader.values())[0])).upper()

  with open("myfile.igd", "rb") as f:
    igd_file = pyigd.IGDReader(f)
    while variant_index < igd_file.num_variants:

      # Collect all the variant indices for the site.
      site_indices = collect_next_site(igd_file, variant_index)
      variant_index = site_indices[-1] + 1

      # Get our position and see if there is an ancestral allele.
      site_position, _, _ = igd_file.get_position_flags_copies(variant_index)
      if site_position < len(ancestral_str):
        flip_ref = ancestral_str[site_position]
        if flip_ref in "ACTG":
          # Get all the alleles
          refs = set(map(lambda i: igd_file.get_ref_allele(i), site_indices))
          alts = list(map(lambda i: igd_file.get_alt_allele(i), site_indices))

          if len(refs) > 1:
              print(f"WARNING: Multiple REF alleles at position {site_position}: {refs}. Skipping site.")
          else:
              old_ref = list(refs)[0]
              if flip_ref == old_ref:
                pass # Nothing to do. REF already matches the ancestral allele.
              elif flip_ref in alts:
                # The old reference samples become a new alternate allele sample list.
                old_ref_samples = get_inverse_sample_list(igd_file, site_indices)
              else:
                pass # Skipped because the ancestral allele was not found in our alternate alleles


IGD Transformation
~~~~~~~~~~~~~~~~~~

:py:meth:`pyigd.IGDWriter` can create arbitrary IGD files, but often we want to just transform one
IGD file into another. Using this type of transformation for filtering is often involved in
bioinformatics pipelines, though that is not the primary use-case target of IGD. Developers of
statistical genetics and population genetics tools often find it useful to manipulate incoming
data for testing and evaluation purposes. For example, adding noise to a dataset to simulate
genotyping error is sometimes used during evaluation, to make the otherwise "clean" simulated
data more realistic. For this reason, we have :py:meth:`pyigd.IGDTransformer`, which makes
simple modification of an existing IGD file very easy.

Filter out low-frequency variants
---------------------------------

Here we transform "file1.igd" into "file2.igd" by removing all of the variants with minor allele
frequency less than ``0.01``. NOTE: this does not filter out the entire site, so if you have multi-
allelic sites there may still be another alternate allele at the site with a frequency exceeding
``0.01``.

::

  import pyigd

  class RemoveLF(pyigd.IGDTransformer):
    def modify_samples(self, position, is_missing, samples, num_copies):
      frequency = len(samples) / self.reader.num_samples
      if frequency < 0.01:
        return None # None means "delete this variant"
      # Otherwise, return the sample list unmodified
      return samples

  with open("file1.igd", "rb") as fin, open("file2.igd", "wb") as fout:
    xformer = RemoveLF(fin, fout, use_bitvectors=False)
    xformer.transform()


Filter out variants outside of range
------------------------------------

Here we transform "file1.igd" into "file2.igd" by only keeping variants within a specified base-
pair range.

::

  import pyigd

  my_range = (5_000_000, 10_000_000)  # 5MBP to 10MBP

  class KeepRange(pyigd.IGDTransformer):
    def modify_samples(self, position, is_missing, samples, num_copies):
      if position >= my_range[0] and position < my_range[1]:
        return samples
      return None # Delete the entire variant

  with open("file1.igd", "rb") as fin, open("file2.igd", "wb") as fout:
    xformer = KeepRange(fin, fout, use_bitvectors=False)
    xformer.transform()


IGD Creation
~~~~~~~~~~~~

IGD creation is done with :py:meth:`pyigd.IGDWriter`.

Writing a genotype matrix
-------------------------

Usually IGD would be useful when you have large datasets, that cannot be loaded into RAM, but
sometimes it is useful to convert a genotype matrix to an IGD file (e.g. for testing, learning
about IGD, etc.). Assume we have a matrix :math:`X` where each of the :math:`N` rows represents
a (haploid) sample and each of the :math:`M` columns represents a bi-allelic variant. Then we
can write this matrix to IGD:

::

  # Matrix is NxM (N = haplotypes, M = mutations)
  def igd_from_matrix(genotype_matrix: numpy.typing.NDArray, filename: str):
    with open(filename, "wb") as fout:
      ploidy = 1
      num_individuals = genotype_matrix.shape[0]

      writer = pyigd.IGDWriter(fout, num_individuals, ploidy)
      # We have to write the header at the start of the file, to "hold its place",
      # even though it doesn't have all the information yet.
      writer.write_header()

      # Write each variant as described in the matrix
      for col in range(genotype_matrix.shape[1]):
        sample_list = numpy.flatnonzero(genotype_matrix[:, col])
        writer.write_variant(col, "0", "1", sample_list)

      # Write the IGD index to the file.
      writer.write_index()

      # Write the variant information (allele strings, mostly)
      writer.write_variant_info()

      # Now we seek back to the beginning of the file and write the header again, since
      # it has been updated with all of the variant information above.
      writer.out.seek(0)
      writer.write_header()

The above is the common pattern for writing an IGD file. Often, people just want to
convert from vcf(.gz) or BGEN, which can be done more efficiently with
`igdtools <https://github.com/aprilweilab/picovcf?tab=readme-ov-file#build-and-run-the-teststools>`_
or `bgen2igd <https://github.com/dcdehaas/bgen2igd>`_, respectively.

IGDWriter can also write identifiers for variants and individuals, via the
:py:meth:`pyigd.IGDWriter.write_variant_ids` and :py:meth:`pyigd.IGDWriter.write_individual_ids`
methods.

Note: the above does a column-wise traversal of the numpy matrix, but it may be more
efficient to first transpose the genotype matrix and then do a row-wise traversal.

Metadata
~~~~~~~~

The metadata is stored separately from the IGD file, in a ``*.meta/`` directory. Each file contains
one metadata array associated with the variants in the IGD file. `numpy.loadtxt <https://numpy.org/doc/2.2/reference/generated/numpy.loadtxt.html>`_
can be used to load the data, and each element in the resulting (1-dimensional) ``numpy.array`` is
the metadata item for that variant index in the IGD file.

For example, consider variant with index ``0 <= j < igd_file.num_variants``. We can use :py:meth:`pyigd.IGDReader.get_position_and_flags`
to get the position of variant ``j``, and given a metadata vector ``meta`` (as loaded by ``numpy.loadtxt``) we can
get that metadata for variant ``j`` via ``meta[j]``.

This is all pretty simple until you modify the IGD file after the metadata has been exported, e.g.
by filtering out some variants. Below is an example that shows you how to match up the variant IDs from
the newly filtered IGD file to the original metadata (which always stores the variant IDs in ``variants.txt``).
This example assumes that you converted a VCF file containing the ``AC`` (allele counts) field from the
``INFO`` column.

::

  import numpy
  import pyigd

  # The original IGD file was called "original.igd", and it was created via: igdtools some.vcf.gz -o original.igd -e all
  original_ids = numpy.loadtxt("original.meta/variants.txt", dtype=str)
  original_ac = numpy.loadtxt("original.meta/info.AC.txt", dtype=int)
  assert len(original_ids) == len(original_ac)

  # Now filter "original.igd" to remove all odd-numbered base-pair positions. This is just a silly example of filtering,
  # to illustrate how metadata is handled after filtering.
  class OddXformer(pyigd.IGDTransformer):
      def modify_samples(self, position, is_missing, samples, num_copies):
          if position % 2 == 1:
              return None
          return samples

  with open("original.igd", "rb") as fin, open("filtered.igd", "wb") as fout:
      xformer = OddXformer(fin, fout)
      xformer.transform()

  # Now open our filtered IGD file, and print out the variant positions plus the corresponding AC value from the metadata
  with open("filtered.igd", "rb") as f:
      igd_file = pyigd.IGDReader(f)
      # The variant IDs that we still have after filtering
      remaining_ids = numpy.array(igd_file.get_variant_ids())

      # Find the indices for the intersection of the variant IDs between the metadata and the filtered IGD.
      _, metadata_indices, remaining_indices = numpy.intersect1d(original_ids, remaining_ids, return_indices=True)
      assert len(metadata_indices) == igd_file.num_variants
      # Sort according to the indices in the filtered IGD file. This way we can just use metadata_indices directly,
      # where metadata_indices[0] is the index into the metadata for the 0th variant in filtered.igd, etc.
      metadata_indices = metadata_indices[numpy.argsort(remaining_indices)]

      # Traverse the filtered variants and their metadata indices together
      for variant_index, metadata_index in zip(range(igd_file.num_variants), metadata_indices):
          assert remaining_ids[variant_index] == original_ids[metadata_index]
          pos, flags = igd_file.get_position_and_flags(variant_index)
          print(f"Position: {pos}, AC={original_ac[metadata_index]}")