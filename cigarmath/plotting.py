"""Various useful plotting operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import numpy as np
import matplotlib.pyplot as plt

from cigarmath.block import reference_block
from cigarmath.block import reference_mapping_blocks


def segments_to_binary(alns, max_genome_size=10_000, deletion_size=50):
    """Given a list of segments create a binary array of all covered regions.
    
    POS0    000000000011111111112222222222333333333
    POS1    012345678901234567890123456789012345678
    REF     AAAAGACCCCCGACTAGCTAGCATGCTATCTAGCTAGCA
    QRY1    AAAAGA-----GAC      
    QRY2                           TGCTA---AGCTAG
    RES     111111000001110000000001111111111111100
    
    alns = [(0, '6M5D3M'),
            (23, '5M3D6M')]
    
    >> segments_to_binary(alns, max_genome_size=38, deletion_size=4)
    """
    
    mapping = np.zeros(max_genome_size, dtype=bool)
    
    for start, cigar in alns:
        block_iter = reference_mapping_blocks(cigar,
                                              reference_start=start,
                                              deletion_split=deletion_size)
        for ref_start, ref_stop in block_iter:
            mapping[ref_start:ref_stop] = 1
        
    return mapping
