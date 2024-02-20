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
    CONSUMES_QUERY,
)





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
