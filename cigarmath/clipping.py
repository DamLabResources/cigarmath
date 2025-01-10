"""Cigar operations involving clipping"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from typing import Tuple, List, Union, TypeVar, Any
from cigarmath.defn import (
    CigarTuples,
    BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP,
    BAM_CMATCH,
    BAM_CEQUAL,
    BAM_CDIFF,
    CONSUMES_QUERY,
    CONSUMES_REFERENCE
)

T = TypeVar('T')  # For generic types in declip function

def left_clipping(cigartuples: CigarTuples, with_hard: bool = True) -> int:
    """Returns the length of clipped bases (hard or soft) on the left side of the alignment

    REF     AAAAACCCCC
    QRY  TTTAAAAACCCCCGGGG
    CGS  SSSMMMMMMMMMMSSSS
    CGT  3S 10M       4S

    >>> left_clipping(cigartuples)
    3
    """
    soft = cigartuples[0][0] == BAM_CSOFT_CLIP
    hard = cigartuples[0][0] == BAM_CHARD_CLIP
    if soft | (with_hard & hard):
        return cigartuples[0][1]
    return 0


def right_clipping(cigartuples: CigarTuples, with_hard: bool = True) -> int:
    """Returns the length of clipped bases (hard or soft) on the right side of the alignment

    REF     AAAAACCCCC
    QRY  TTTAAAAACCCCCGGGG
    CGS  SSSMMMMMMMMMMSSSS
    CGT  3S 10M       4S

    >>> right_clipping(cigartuples)
    4
    """
    soft = cigartuples[-1][0] == BAM_CSOFT_CLIP
    hard = cigartuples[-1][0] == BAM_CHARD_CLIP
    if soft | (with_hard & hard):
        return cigartuples[-1][1]
    return 0


def declip(cigartuples: CigarTuples, *args: T) -> CigarTuples:
    """Return a set of cigartuples with clipping removed, if any.

    REF     AAAAACCCCC
    QRY  TTTAAAAACCCCCGGGG
    CGS  SSSMMMMMMMMMMSSSS
    CGT  3S 10M       4S

    >>>> declip(cigartuples)
    (0, 10)
    
    You can also provide sequences or lists as extra *args
    which will be clipped based on the tuples.
    
    >>>> seq = 'TTTAAAAACCCCCGGGG'
    >>>> tups, clipped_seq = declip(cigartuples, seq)
    (0, 10)
    AAAAACCCCC
    """
    left_clip = (cigartuples[0][0] == BAM_CSOFT_CLIP) | (
        cigartuples[0][0] == BAM_CHARD_CLIP
    )
    right_clip = (cigartuples[-1][0] == BAM_CSOFT_CLIP) | (
        cigartuples[-1][0] == BAM_CHARD_CLIP
    )
    
    cigarstart = 1 if left_clip else 0
    cigarend = -1 if right_clip else None
    
    if args:
        clipstart = cigartuples[0][1] if left_clip else 0
        clipend = -cigartuples[-1][1] if right_clip else None
        print(clipstart, clipend)
        clipped_args = [items[clipstart:clipend] for items in args]
        return cigartuples[cigarstart:cigarend], *clipped_args
    
    return cigartuples[cigarstart:cigarend]


def is_hard_clipped(cigartuples: CigarTuples) -> bool:
    """Return True if the cigar indicates HARD clipping

    REF     AAAAACCCCC
    QRY     AAAAACCCCC
    CGS  HHHMMMMMMMMMMHHHH
    CGT  3H 10M       4H

    >>>> is_hard_clipped(cigartuples)
    True
    """
    return (cigartuples[0][0] == BAM_CHARD_CLIP) or (
        cigartuples[-1][0] == BAM_CHARD_CLIP
    )


def softclipify(cigartuples: CigarTuples, required_mapping: int = 1) -> Tuple[CigarTuples, int]:
    """Converts initial and final cigars into softclips
    
    REF    --AAAAGACCCCCGACTCGTTA---
    QUE    TT----AACCCCCGAC----TAGCA
    CIG    IIDDDDMMMMMMMMMMDDDDMMIII
        
    OUT    SS    MMMMMMMMMMDDDDMMSSS required_mapping = 1
    OUT    SS    MMMMMMMMMM    SSSSS required_mapping = 4
    """
    left_ind, left_soft_sz = _decide_softclip_end(cigartuples, required_mapping)
    right_ind, right_soft_sz = _decide_softclip_end(cigartuples[::-1], required_mapping)

    offset = 0
    if left_soft_sz:
        for op, sz in cigartuples[:left_ind]:
            if op in CONSUMES_REFERENCE:
                offset += sz
    
    left = [(BAM_CSOFT_CLIP, left_soft_sz)] if left_soft_sz else []
    right = [(BAM_CSOFT_CLIP, right_soft_sz)] if right_soft_sz else []
    
    middle_slc = slice(left_ind, -right_ind or None)
    
    return left+cigartuples[middle_slc]+right, offset


def _decide_softclip_end(cigartuples: CigarTuples, required_mapping: int) -> Tuple[Union[int, None], Union[int, None]]:
    """Helper function to search cigartuples until you find a MAPPING block of a particular size.
    
    REF    --AAAAGACCCCCGACTCGTTA---
    QUE    TT----AACCCCCGAC----TAGCA
    CIG    IIDDDDMMMMMMMMMMDDDDMMIII 
    
    cigartuples = [(BAM_CINS, 2), (BAM_CDEL, 4), (BAM_CMATCH, 10), 
                   (BAM_CDEL, 4), (BAM_CMATCH, 2), (BAM_CINS, 3)]
    """
    MAPPING_OP = set([BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF])
    
    soft_sz = 0
    
    for num, (op, sz) in enumerate(cigartuples):
        if (op in MAPPING_OP) and (sz >= required_mapping):
            return num, soft_sz
        elif op in CONSUMES_QUERY:
            soft_sz += sz

    return None, None



# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
