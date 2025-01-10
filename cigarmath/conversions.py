"""Ways to convert between alignment formats"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from typing import List, Tuple, Iterator, Optional, Union
from cigarmath.block import reference_block
from cigarmath.block import reference_mapping_blocks
from cigarmath.clipping import left_clipping
from cigarmath.clipping import right_clipping
from cigarmath.clipping import declip
from cigarmath.clipping import softclipify
from cigarmath.cigarmath import collapse_adjacent_blocks
from cigarmath.defn import (
    CigarTuples,
    CONSUMES_REFERENCE,
    CONSUMES_QUERY,
    BAM_CMATCH,
    BAM_CINS,
    BAM_CDEL,
    BAM_CEQUAL,
    BAM_CDIFF,
    BAM_CSOFT_CLIP,
)

# Type aliases
AlignmentSegment = Tuple[int, str]
AlignmentPair = Tuple[int, Optional[int]]
MSA = str  # Multiple Sequence Alignment string

def segments_to_binary(
    alns: List[AlignmentSegment], 
    max_genome_size: int = 10_000, 
    deletion_size: int = 50, 
    mapping: Optional[List[bool]] = None
) -> List[bool]:
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
        mapping = [False] * max_genome_size
    
    for start, cigar in alns:
        block_iter = reference_mapping_blocks(cigar,
                                           reference_start=start,
                                           deletion_split=deletion_size)
        for ref_start, ref_stop in block_iter:
            mapping[ref_start:ref_stop] = [True] * (ref_stop-ref_start)
        
    return mapping


def cigartuples2pairs(
    cigartuples: CigarTuples, 
    reference_start: int = 0, 
    verbose: bool = False, 
    clipping_fill: Optional[int] = None
) -> List[AlignmentPair]:
    
    reference: List[Optional[int]] = []
    query: List[int] = []
    
    left_clip = left_clipping(cigartuples)
    if left_clip:
        if verbose: print('Adding left clipping', left_clip)
        query += list(range(left_clip+1))
        reference = [clipping_fill] * left_clip
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


def msa2cigartuples(ref_msa: MSA, query_msa: MSA) -> Tuple[int, CigarTuples]:
    """Given a pair of multiple-sequence alignments return reference_start and cigartuples.
    
    REF     AAAAGACCCCCGACTAGCTAGCATGCT----ATCTAGCTAGCA
    QRY     ----AACCCCCGAC----TAGCATGCTTTTTATCTAGCT----
    CIGAR       MMMMMMMMMMDDDDMMMMMMMMIIIIIMMMMMMMM
    
    >> ref_start, cigartuples = msa2cigartuples(ref_msa, query_msa)
    ref_start = 4
    cigartuples = [(0, 10), (2, 4), (0, 9), (1, 4), (0, 8)]
    """
    
    assert len(ref_msa) == len(query_msa)
    
    cigartuples: CigarTuples = []
    cur_op: Optional[int] = None
    cur_sz = 0
    reference_start: Optional[int] = None
    offset = 0
    query_started = False
    
    non_double = ((r, q) for r, q in zip(ref_msa, query_msa) if (r!='-') | (q!='-'))
    
    for n, (r, q) in enumerate(non_double):
        query_started |= q != '-'
        # Don't start until the first non-gap query character
        if query_started:
            this_op = _decide_op(r, q)
            # print(n, r, q, this_op)
            if (cur_op is None) and (reference_start is None):
                # First OP
                cur_op, cur_sz = this_op, 1
                reference_start = n
            elif this_op == cur_op:
                # Extend this op
                cur_sz += 1
            else:
                # New op
                cigartuples.append((cur_op, cur_sz))
                cur_op, cur_sz = this_op, 1
        else:
            offset += 1
    
    if cur_op is not None:
        cigartuples.append((cur_op, cur_sz))
    
    collapsed_tuples = collapse_adjacent_blocks(cigartuples)
    softclipped_tuples, clip_offset = softclipify(cigartuples, required_mapping=1)

    assert reference_start is not None
    return reference_start + clip_offset, softclipped_tuples
    
    


    
def cigartuples2msa(
    reference: str, 
    query: str, 
    reference_start: int, 
    cigartuples: CigarTuples
) -> Tuple[MSA, MSA]:
    """Converts a cigartuple and reference start into a Multiple Sequence Alignment"""
    raise NotImplementedError
    
    
    
    
    
    
def _decide_op(ref_letter: str, query_letter: str, extended: bool = False) -> int:
    """Determines the cigarop from reference and query MSA positions"""
    
    if (ref_letter == '-') & (query_letter == '-'):
        # Both are gaps
        return None
    elif ref_letter == '-':
        # Reference gap, therefore insertion
        return BAM_CINS
    elif query_letter == '-':
        # query gap, therefore deletion
        return BAM_CDEL
    elif ref_letter == query_letter:
        # identical
        return BAM_CEQUAL if extended else BAM_CMATCH
    elif ref_letter != query_letter:
        # mismatch
        return BAM_CDIFF if extended else BAM_CMATCH
    
    raise AssertionError('Did not expect to get here') 