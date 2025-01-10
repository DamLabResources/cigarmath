"""Dead-simple io"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import random
from typing import Union, Iterator, Optional, Tuple, TYPE_CHECKING
from cigarmath.defn import CigarTuples

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