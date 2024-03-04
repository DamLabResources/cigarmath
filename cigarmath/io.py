"""Dead-simple io"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import random


def segment_stream_pysam(path, mode='rt', fetch=None, 
                         min_mapq=0, downsample = None,
                         as_tuples=False):
    """
    Yield AlignedSegments from sam/bam file with pysam.
    """
    
    import pysam
    
    # Open the SAM/BAM file with PySAM
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
                    
                    
                    
def _downsample(stream, frac):
    
    for item in stream:
        if random.random() < frac:
            yield item