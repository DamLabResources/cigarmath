"""Various useful cigar operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import CONSUMES_QUERY



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


def inferred_reference_length(cigartuples):
    """Returns the expected length of query_sequence based on cigartuples
    
    POS  01234567890  12345
    REF     AAAAGACC--CCC
    QRY     AAAA-ACCGGCCC
    CGS  HHHMMMMDMMMIIMMMHHHH
    CGT  3H 4M 1D3M 2I 3M 4H
    
    >>>> inferred_reference_length(cigartuples)
    11
    """

    return sum(block_size for bam_num, block_size in cigartuples if bam_num in CONSUMES_REFERENCE)