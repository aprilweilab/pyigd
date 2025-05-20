.. pyigd documentation master file, created by
   sphinx-quickstart on Wed Aug 21 14:20:22 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyigd Documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   IGD Overview <igd_overview>
   Example Usage <examples>
   API Reference <pyigd>
   igdtools

pyigd is a Python library for reading and writing Indexable Genotype Data
(IGD) files. `IGD is a binary format <https://github.com/aprilweilab/picovcf/blob/main/IGD.FORMAT.md>`_ that is very simple and uses no
compression, but tends to be quite small and fast.

IGDReader example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example of how to read an IGD file. See :py:class:`pyigd.IGDReader` for the relevant API reference.

::

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
          print(f"REF: {igd_file.get_ref_allele(variant_index)}, ALT: {igd_file.get_alt_allele(variant_index)}")
          position, is_missing, sample_list = igd_file.get_samples(variant_index)
          print( (position, is_missing, len(sample_list))  )

IGDWriter example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example of how to write an IGD file. See :py:class:`pyigd.IGDWriter` for the relevant API reference. Some of the methods invoked have a required order, and in general all the variants (sample info) needs to be written immediately after the header. The IGDWriter maintains some state, and for large datasets can use a non-trivial amount of RAM to maintain indexes while writing variant data.

::

  num_individuals = 100
  with open("output.igd", "wb") as f:
      w = IGDWriter(f, num_individuals)
      w.write_header()
      w.write_variant(100, "A", "G", [1, 22, 99, 101])
      w.write_variant(101, "A", "C", list(range(200)))
      w.write_index()
      w.write_variant_info()
      w.write_individual_ids([]) # optional
      w.write_variant_ids([])    # optional
      # We write the header twice because the IGDWriter is updating fields
      # on the header as we call the other methods above. For example, it is
      # counting the number of variants being written and that gets updated
      # in the header the second time we write it.
      f.seek(0)
      w.write_header()


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
