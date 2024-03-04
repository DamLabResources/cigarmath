"""Ways to convert between alignment formats"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"



from cigarmath.block import reference_block
from cigarmath.block import reference_mapping_blocks

from cigarmath.clipping import left_clipping
from cigarmath.clipping import right_clipping
from cigarmath.clipping import declip

from cigarmath.defn import CONSUMES_REFERENCE
from cigarmath.defn import CONSUMES_QUERY


def segments_to_binary(alns, max_genome_size=10_000, deletion_size=50, mapping=None):
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
    
    if mapping is None:
        mapping = [False]*max_genome_size
    
    for start, cigar in alns:
        block_iter = reference_mapping_blocks(cigar,
                                              reference_start=start,
                                              deletion_split=deletion_size)
        for ref_start, ref_stop in block_iter:
            mapping[ref_start:ref_stop] = [True]*(ref_stop-ref_start)
        
    return mapping


def cigartuples2pairs(cigartuples, reference_start = 0, verbose=False, clipping_fill = None):
    
    reference = []
    query = []
    
    left_clip = left_clipping(cigartuples)
    if left_clip:
        if verbose: print('Adding left clipping', left_clip)
        query += list(range(left_clip+1))
        reference = [clipping_fill]*left_clip
        if verbose: print('Currently', list(zip(query, reference)))
    
    for op, sz in declip(cigartuples):
        qstart = query[-1] if query else -1
        if verbose: print('Considering', op, sz)
        if op in CONSUMES_QUERY:
            if verbose: print('Consumes query, adding', list(range(qstart+1, qstart+sz+1)))
            query += list(range(qstart+1, qstart+sz+1))
        else:
            if verbose: print('Skips query, adding', [qstart]*sz)
            query += [query[-1]]*sz
        
        rstart = reference[-1] if (reference and reference[-1]) else reference_start-1
        if op in CONSUMES_REFERENCE:
            if verbose: print('Consumes Ref, adding', list(range(rstart+1, rstart+sz+1)))
            
            reference += list(range(rstart+1, rstart+sz+1))
        else:
            if verbose: print('Skips Ref, adding',  [rstart]*sz)
            reference += [rstart]*sz
        
        if verbose: print('Currently', list(zip(query, reference)))
            
    right_clip = right_clipping(cigartuples)
    if right_clip:
        if verbose: print('Adding right clipping:', right_clip)
        query += list(range(query[-1]+1, query[-1]+right_clip+1))
        reference += [clipping_fill]*right_clip
        if verbose: print('Currently', list(zip(query, reference)))
            
    return list(zip(query, reference))