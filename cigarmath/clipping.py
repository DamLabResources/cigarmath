"""Cigar operations involving clipping"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import BAM_CSOFT_CLIP
from cigarmath.defn import BAM_CHARD_CLIP
from cigarmath.defn import BAM_CMATCH
from cigarmath.defn import BAM_CEQUAL
from cigarmath.defn import BAM_CDIFF
from cigarmath.defn import CONSUMES_QUERY
from cigarmath.defn import CONSUMES_REFERENCE


def left_clipping(cigartuples, with_hard=True):
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


def right_clipping(cigartuples, with_hard=True):
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


def declip(cigartuples, *args):
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
    
    cigarstart = 1 if left_clip else 0  # Start from the second element if there's a left clip, otherwise start from the first
    cigarend = -1 if right_clip else None  # End one element early if there's a right clip, otherwise go to the end
    
    if args:
        clipstart = cigartuples[0][1] if left_clip else 0
        clipend = -cigartuples[-1][1] if right_clip else None
        print(clipstart, clipend)
        clipped_args = [items[clipstart:clipend] for items in args]
        return cigartuples[cigarstart:cigarend], *clipped_args
    
    return cigartuples[cigarstart:cigarend]
        
    
    

    


def is_hard_clipped(cigartuples):
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


def softclipify(cigartuples, required_mapping = 1):
    """Converts initial and final cigars into softclips
    
    REF    --AAAAGACCCCCGACTCGTTA---
    QUE    TT----AACCCCCGAC----TAGCA
    CIG    IIDDDDMMMMMMMMMMDDDDMMIII
        
    OUT    SS    MMMMMMMMMMDDDDMMSSS required_mapping = 1
    OUT    SS    MMMMMMMMMM    SSSSS required_mapping = 4
    
    
    This is useful when converting MSA into cigartuples where there may be leading 'insertions'
    before the initial reference match leading to invalid cigarstrings.
    
    """
    
    left_ind, left_soft_sz = _decide_softclip_end(cigartuples, required_mapping)
    right_ind, right_soft_sz = _decide_softclip_end(cigartuples[::-1], required_mapping)

    # The new offset is the sum of reference consuming blocks
    # that are being discarded
    offset = 0
    if left_soft_sz:
        for op, sz in cigartuples[:left_ind]:
            if op in CONSUMES_REFERENCE:
                offset += sz
    
    # Create new softclip blocks
    left = [(BAM_CSOFT_CLIP, left_soft_sz)] if left_soft_sz else []
    right = [(BAM_CSOFT_CLIP, right_soft_sz)] if right_soft_sz else []
    
    # Get a slice for the rematining blocks
    middle_slc = slice(left_ind, -right_ind or None)
    
    return left+cigartuples[middle_slc]+right, offset


MAPPING_OP = set([BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF])

def _decide_softclip_end(cigartuples, required_mapping):
    """Helper function to search cigartuples until you find a MAPPING block of a particular size.
    
    REF    --AAAAGACCCCCGACTCGTTA---
    QUE    TT----AACCCCCGAC----TAGCA
    CIG    IIDDDDMMMMMMMMMMDDDDMMIII 
    
    cigartuples = [(BAM_CINS, 2), (BAM_CDEL, 4), (BAM_CMATCH, 10), 
                   (BAM_CDEL, 4), (BAM_CMATCH, 2), (BAM_CINS, 3)]
    
    ind, soft_sz = _decide_softclip_end(cigartuples, 1)
    ind = 2
    soft_sz = 2
    
    """
    
    soft_sz = 0
    
    for num, (op, sz) in enumerate(cigartuples):
        if (op in MAPPING_OP) and (sz >=required_mapping):
            return num, soft_sz
        elif op in CONSUMES_QUERY:
            soft_sz += sz

    return None, None



# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
