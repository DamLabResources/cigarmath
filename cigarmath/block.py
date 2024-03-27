"""Various useful block-level operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import CONSUMES_REFERENCE
from cigarmath.defn import CONSUMES_QUERY
from cigarmath.defn import BAM_CSOFT_CLIP
from cigarmath.defn import BAM_CHARD_CLIP
from cigarmath.defn import BAM_CDEL
from cigarmath.defn import BAM_CREF_SKIP
from cigarmath.clipping import left_clipping


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
    return sum(
        block_size
        for bam_num, block_size in cigartuples
        if bam_num in CONSUMES_REFERENCE
    )


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


def query_offset(cigartuples):
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
    consumes_query_offset.discard(BAM_CSOFT_CLIP), consumes_query_offset.discard(
        BAM_CHARD_CLIP
    )

    # add up remaining query-consuming blocks
    return sum(
        block_size
        for bam_num, block_size in cigartuples
        if bam_num in consumes_query_offset
    )


def query_start(cigartuples):
    """Return the start position on the query of this alignment.

    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H

    >>>> query_start(cigartuples)
    3
    """
    
    return left_clipping(cigartuples)


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
    offset = query_offset(cigartuples)

    return query_start, query_start + offset


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


def reference_mapping_blocks(cigartuples, reference_start=0, deletion_split=10):
    """Yield (reference_start, reference_stop) blocks of mapped sites split by deletions.
    
    POS0  000000000011111111112222222222
    POS1  012345678901234567890123456789

    CGS      MMMMMMMDDDMMMMDDDDDDMMMM
    
    >> reference_mapping_blocks(cigartuples, reference_start=3, deletion_split=5)
    (3, 16)
    (22, 26)
    """
    
    del_ops = {BAM_CDEL, BAM_CREF_SKIP}
    
    left, right = reference_start, reference_start
    for op, sz in cigartuples:
        if (op in del_ops) and (sz >= deletion_split):
            yield (left, right)
            left = right + sz
        right += sz * (op in CONSUMES_REFERENCE)
        
    yield left, right


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


# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
