"""Cigar operations involving clipping"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import BAM_CSOFT_CLIP
from cigarmath.defn import BAM_CHARD_CLIP

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


def declip(cigartuples):
    """Return a set of cigartuples with clipping removed, if any.
    
    REF     AAAAACCCCC
    QRY  TTTAAAAACCCCCGGGG
    CGS  SSSMMMMMMMMMMSSSS
    CGT  3S 10M       4S
    
    >>>> declip(cigartuples)
    (0, 10)
    """
    
    left_clip = (cigartuples[0][0] == BAM_CSOFT_CLIP) | (
        cigartuples[0][0] == BAM_CHARD_CLIP
    )
    right_clip = (cigartuples[-1][0] == BAM_CSOFT_CLIP) | (
        cigartuples[-1][0] == BAM_CHARD_CLIP
    )

    if left_clip and right_clip:
        return cigartuples[1:-1]
    elif right_clip:
        return cigartuples[:-1]
    elif left_clip:
        return cigartuples[1:]
    return cigartuples


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


# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved