"""
Module for reading and writing Indexable Genotype Data (IGD) files.
"""

from .readwrite import (  # noqa: F401
    BpPosFlags,
    IGDReader,
    IGDConstants,
    IGDHeader,
    IGDTransformer,
    IGDWriter,
    flags_is_missing,
)
