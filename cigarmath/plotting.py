"""Various useful plotting operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import numpy as np
import matplotlib.pyplot as plt

from cigarmath.block import reference_block
from cigarmath.block import reference_deletion_blocks


def segments_to_binary(cigartuples, reference_starts, max_genome_size=10_000, deletion_size=50):
    """Given a list of segments create a binary array of all covered regions.
    
    REF     AAAAGACCCCCGACTAGCTAGCATGCTATCTAGCTAGCA
    QRY1    AAAAGA-----GAC      
    QRY2                           TGCTA---AGCTAG
    RES     111111000001110000000001111100011111100
    
    """
    
    mapping = np.zeros(max_genome_size, dtype=bool)
    
    for cigar, start in zip(cigartuples, reference_starts):
        
        ref_start, ref_stop = reference_block(cigar, start)
        mapping[ref_start:ref_stop] = 1
        
        del_blocks = reference_deletion_blocks(cigar, 
                                               reference_start = start,
                                               min_size = deletion_size)
        for del_begin, del_end in del_blocks:
            mapping[del_begin:del_end] = 0
            
    return mapping
        