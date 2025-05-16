"""
Extra functionality that is not core to accessing an IGD, but helpful for manipulating
the information in/related to an IGD.
"""

from typing import List

try:
    import numpy
except ImportError:
    numpy = None  # type: ignore


def collect_next_site(igd_reader, variant_index: int) -> List[int]:
    """
    Given a variant index to start at, iterate all consecutive variants that have the same position
    and return the variant indices for that position. The returned indices are in ascending order.

    :param igd_reader: The IGDReader representing the IGD file.
    :type igd_reader: IGDReader
    :param variant_index: The variant index to start scanning from.
    :type variant_index: int
    :return: The list of subsequent variant indexes that all share the same position.
    :rtype: List[int]
    """
    results = []
    site_position = None
    while variant_index < igd_reader.num_variants:
        position, _, _ = igd_reader.get_position_flags_copies(variant_index)
        if site_position is not None and position != site_position:
            break
        site_position = position
        results.append(variant_index)
        variant_index += 1
    return results


# For unphased data, we can have an individual which has 1 copy of A, 2 copies of G, and a ploidy of 4.
# Then the reference has 1 copy (the inverse). We need to look @ all copies to determine ploidy copies
def get_inverse_sample_list(igd_reader, variant_indices: List[int]) -> List[List[int]]:
    """
    Given a list of variant indices, compute the set of samples that are covered by those indices, and then
    invert that list. Works for both phased and unphased data. One usage for this function is to get an explicit
    list of samples that have the reference, since the reference is stored implicitly in an IGD.

    :param igd_reader: The IGDReader representing the IGD file.
    :type igd_reader: IGDReader
    :param variant_indices: The indices of the variants who total sample set (i.e., all of them unioned together)
        you want the inverse of.
    :type variant_indices: List[int]
    :return: A list of sample lists, one for each number of copies between 1...PLOIDY. For phased data, there
        will always be a single sample list. For unphased data, there will always be PLOIDY sample lists.
    :rtype: List[List[int]]
    """
    # If you are inverting a sample list, then sparseness no longer matters. Either:
    # - The original list was not sparse, so you have to iterate a lot of items.
    # - The original list was sparse, so the inverse is not.
    # Hence, the inversion is a non-sparse operation and we use a matrix here.
    assert numpy is not None, "numpy required for this operation: 'pip install numpy'"
    ploidy = igd_reader.ploidy
    coverage_count = numpy.full(
        igd_reader.num_samples, 1 if igd_reader.is_phased else ploidy, dtype=numpy.int32
    )
    for i in variant_indices:
        _, _, num_copies = igd_reader.get_position_flags_copies(i)
        # For phased data.
        if num_copies == 0:
            num_copies = 1
        _, _, samples = igd_reader.get_samples(i)
        this_vector = numpy.zeros(igd_reader.num_samples, dtype=numpy.int32)
        numpy.put(this_vector, samples, num_copies)
        coverage_count -= this_vector
    result = []
    # If we have negative values, it means that the input data had overlapping variants that added up to
    # more than PLOIDY for some samples. In real data this often happens at sites containing indels/SVs,
    # where the overlap isn't very well resolved and so a particular sample is called as having SNV "x" and
    # indel "y", which isn't always very interpretable.
    # We return no results for any such sites.
    if numpy.all(coverage_count >= 0):
        if igd_reader.is_phased:
            result.append(numpy.flatnonzero(coverage_count).tolist())
        else:
            for num_copies in range(1, ploidy + 1):
                result.append(numpy.flatnonzero(coverage_count == num_copies).tolist())
    return result
