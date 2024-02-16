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
    cigartuples2cigarstring
)

CONSUMES_REFERENCE = {bam_num for bam_num, cigar_letter in enumerate(NTS) if cigar_letter.consumes_ref}
CONSUMES_QUERY = {bam_num for bam_num, cigar_letter in enumerate(NTS) if cigar_letter.consumes_query}



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


def reference_offset(cigartuples):
    """Calculate the length of the reference mapping block based on cigartuples.
    
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> reference_offset(cigartuples)
    10
    """

    # add up reference-consuming blocks
    return sum(block_size for bam_num, block_size in cigartuples if bam_num in CONSUMES_REFERENCE)


def reference_block(cigartuples, reference_start=0):
    """Returns a tuple of the reference (start, end) positions of the aligned segment
    
    POS  01234567890  12345
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> reference_block(cigartuples, reference_start=3)
    (3, 13)
    """

    offset = reference_offset(cigartuples)
    return reference_start, reference_start + offset


def query_offest(cigartuples):
    """Calculate the length of the query mapping block based on cigartuples.
    
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> query_offset(cigartuples)
    11
    """
    
    consumes_query_offset = set(CONSUMES_QUERY)
    # remove clipping
    consumes_query_offset.discard(BAM_CSOFT_CLIP), consumes_query_offset.discard(BAM_CHARD_CLIP)

    # add up remaining query-consuming blocks
    return sum(block_size for bam_num, block_size in cigartuples if bam_num in consumes_query_offset)


def query_block(cigartuples):
    """Returns a tuple of the query (start, end) positions of the aligned segment
    
    POS  01234567890  12345
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> query_block(cigartuples)
    (3, 12)
    """

    # check clipping
    query_start = left_clipping(cigartuples)

    # add up remaining query-consuming blocks
    offset = query_offest(cigartuples)

    return query_start, query_start + offset


def inferred_query_sequence_length(cigartuples):
    """Returns the expected length of query_sequence based on cigartuples
    
    POS  01234567890  12345
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> inferred_query_sequence_length(cigartuples)
    19
    """

    return sum(block_size for bam_num, block_size in cigartuples if bam_num in CONSUMES_QUERY)


def block_overlap_length(block_a, block_b):
    """Given two (start,stop) tuples, return their overlap.
    negative values indicate distance to overlap.
    
    POS0  000000000011111111112222222222
    POS1  012345678901234567890123456789
    
    BLKA      ----------
    BLKB            ----------
    
    >>> block_overlap_length((4, 13), (10, 19))
    4
    
    POS0  000000000011111111112222222222
    POS1  012345678901234567890123456789
    
    BLKA      -----
    BLKB                  ----------
    
    >>> block_overlap_length((4, 8), (16, 25))
    -7
    """

    return min(block_a[1], block_b[1]) - max(block_a[0], block_b[0])


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