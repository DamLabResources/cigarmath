"""Top-level package for Cigar Math."""

__author__ = """Will Dampier"""
__email__ = "wnd22@drexel.edu"
__version__ = "0.1.0"


from .clipping import left_clipping
from .clipping import right_clipping
from .clipping import declip
from .clipping import is_hard_clipped
from .clipping import left_clipping
from .clipping import softclipify


from .block import reference_offset
from .block import reference_block
from .block import query_start
from .block import query_offset
from .block import query_block
from .block import block_overlap_length
from .block import reference_offset
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
from .conversions import softclipify

from .cigarmath import collapse_adjacent_blocks

from .mapping import cigar_iterator
from .mapping import reference2query
from .mapping import query2reference