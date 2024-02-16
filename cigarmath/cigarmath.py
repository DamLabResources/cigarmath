"""Various useful cigar operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import (
    BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP,
    NTS,
    BAM_CDEL,
    BAM_CREF_SKIP,
    CIGAR_HDRS,
    CIGAR2BAM,
    CONSUMES_REFERENCE,
    CONSUMES_QUERY
)
















def reference_deletion_blocks(cigartuples, reference_start=0, min_size=1):
    """Yield (reference_start, reference_stop) blocks of deletions larger than minimum size
    
    POS0  000000000011111111112222222222
    POS1  012345678901234567890123456789
    
    CGS      MMMMDDDDDDMMMMDDDDDDMMMM 
    
    >>> reference_deletion_blocks(cigartuples, reference_start=3)
    (7, 12)
    (17, 22)
    """

    del_ops = {BAM_CDEL, BAM_CREF_SKIP}
    for op, sz in cigartuples:
        if (op in del_ops) and (sz >= min_size):
            yield reference_start, reference_start + sz
        reference_start += (op in CONSUMES_REFERENCE) * sz


def simplify_blocks(cigartuples, collapse=True):
    """Replace extended cigars (= and X) with M.
    
    POS  0123456789012345
    REF     AAAAGACCCCC
    QRY     AAAAACCGGCC
    CGS     ====xx=xx==
    
    >>>> simplify_blocks(cigartuples)
    [(0, 11)]
    """

    if collapse:
        simpled = simplify_blocks(cigartuples, collapse=False)
        yield from collapse_adjacent_blocks(simpled)
    else:
        for op, sz in cigartuples:
            if (op == 7) or (op == 8):
                op = 0
            yield op, sz


def collapse_adjacent_blocks(cigartuples):
    """Merge identical adjacent blocks.
    For example:
    3M3M2I -> 6M2I
    """

    last_op = None
    last_sz = 0
    for op, sz in cigartuples:
        if last_op is None:
            # first block
            last_op = op
            last_sz = sz

        elif last_op == op:
            # matching block
            last_sz += sz

        else:
            # new block

            # yield the previous
            yield last_op, last_sz

            # start a new one
            last_op, last_sz = op, sz

    # at the end of the loop, yield what's left
    yield last_op, last_sz

    
# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved