"""Dead-simple io"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import random
from itertools import groupby
from typing import Union, Iterator, Optional, Tuple, TYPE_CHECKING, List
from cigarmath.defn import CigarTuples
from cigarmath.combine import combine_multiple_alignments

if TYPE_CHECKING:
    try:
        import pysam
    except ImportError:
        pass

def segment_stream_pysam(
    path: str,
    mode: str = 'rt',
    fetch: Optional[str] = None,
    min_mapq: int = 0,
    downsample: Optional[float] = None,
    as_tuples: bool = False
) -> Iterator[Union[Tuple[int, CigarTuples], "pysam.AlignedSegment"]]:
    """
    Yield AlignedSegments from sam/bam file with pysam.
    """

    import pysam

    with pysam.AlignmentFile(path, mode, check_sq=False) as samfile:
        iterable = samfile
        
        if fetch:
            iterable = iterable.fetch(fetch)
            
        if downsample:
            iterable = _downsample(iterable, downsample)
        
        for segment in iterable:
            if segment.mapping_quality > min_mapq:
                if as_tuples:
                    if segment.cigartuples:
                        yield (segment.reference_start, segment.cigartuples)
                else:
                    yield segment


def _downsample(stream: Iterator, frac: float) -> Iterator:
    """Randomly sample items from a stream with given fraction."""
    for item in stream:
        if random.random() < frac:
            yield item


def _combine_aligned_segments(segments: Iterator) ->Tuple[int, CigarTuples, List]:
    """Combine aligned segments into a single alignment."""
    
    segments = list(segments)
    try:
        packed = [(segment.reference_start, segment.cigartuples) for segment in segments]
        new_start, new_cigars = combine_multiple_alignments(packed)
    except ValueError:
        return None, None, segments
    return new_start, new_cigars, segments


def _get_primary_segment(segments: List):
    """Get the primary segment from a list of aligned segments."""
    for segment in segments:
        if segment.is_primary_alignment:
            return segment
    return segments[0]

def combined_segment_stream(segments: Iterator) -> Iterator[Tuple[int, CigarTuples, List]]:
    """Combine aligned segments into a single alignment."""

    for _, segments in groupby(segments, key=lambda x: x.query_name):
        segments = list(segments)
        if len(segments) > 1:
            try:
                new_start, new_cigars, segments = _combine_aligned_segments(segments)
                yield (new_start, new_cigars, segments)
            except ValueError:
                primary = _get_primary_segment(segments)
                yield (primary.reference_start, primary.cigartuples, segments)
        else:
            yield (segments[0].reference_start, segments[0].cigartuples, segments)

    