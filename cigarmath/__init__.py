"""Top-level package for Cigar Math."""

import importlib.metadata
import importlib.util

__author__ = """Will Dampier"""
__email__ = "wnd22@drexel.edu"

from .clipping import left_clipping
from .clipping import right_clipping
from .clipping import declip
from .clipping import is_hard_clipped
from .clipping import softclipify


from .block import reference_offset
from .block import reference_block
from .block import query_start
from .block import query_offset
from .block import query_block
from .block import block_overlap_length
from .block import reference_mapping_blocks
from .block import reference_deletion_blocks

from .inference import inferred_query_sequence_length
from .inference import inferred_reference_length

from .defn import cigarstr2tup
from .defn import cigartup2str

from . import io

from .conversions import segments_to_binary
from .conversions import cigartuples2pairs

from .conversions import msa2cigartuples

from .cigarmath import collapse_adjacent_blocks

from .mapping import reference2query
from .mapping import query2reference

from .iterators import cigar_iterator
from .iterators import cigar_iterator_reference_slice
from .iterators import liftover
from .iterators import iterator_attach

from .combine import combine_multiple_alignments
from .combine import combine_adjacent_alignments
from .combine import trim_alignment

from .pileup import depth

from .rearrangement import infer_rearrangements
from .rearrangement import format_read_rearrangement_summary
from .rearrangement import rearrangement_segment_stream
from .rearrangement import reference_lengths_from_pysam_header
from .rearrangement import RearrangementEvent


def _read_runtime_version() -> str:
    try:
        return importlib.metadata.version("cigarmath")
    except importlib.metadata.PackageNotFoundError:
        if importlib.util.find_spec("cigarmath") is not None:
            return "0.0.0+uninstalled"
        raise


__version__ = _read_runtime_version()

__all__ = [
    "left_clipping",
    "right_clipping",
    "declip",
    "is_hard_clipped",
    "softclipify",
    "reference_offset",
    "reference_block",
    "query_start",
    "query_offset",
    "query_block",
    "block_overlap_length",
    "reference_mapping_blocks",
    "reference_deletion_blocks",
    "inferred_query_sequence_length",
    "inferred_reference_length",
    "cigarstr2tup",
    "cigartup2str",
    "io",
    "segments_to_binary",
    "cigartuples2pairs",
    "msa2cigartuples",
    "collapse_adjacent_blocks",
    "reference2query",
    "query2reference",
    "cigar_iterator",
    "cigar_iterator_reference_slice",
    "liftover",
    "iterator_attach",
    "combine_multiple_alignments",
    "combine_adjacent_alignments",
    "trim_alignment",
    "depth",
    "infer_rearrangements",
    "format_read_rearrangement_summary",
    "rearrangement_segment_stream",
    "reference_lengths_from_pysam_header",
    "RearrangementEvent",
    "__version__",
]
