"""Create mappings between read & reference based on CIGARs"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from typing import Iterator, Optional, Tuple
from cigarmath.defn import CigarTuples
from cigarmath.iterators import cigar_iterator

def reference2query(cigartuples: CigarTuples, reference_start: int = 0) -> Iterator[Optional[int]]:
    """Create a generator the same size as the reference alignment
    that maps positions in the reference to positions in the query.
    
    RPOS    0123  456789 # Index within the reference
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query 
    
    r2q =     12  5NN6  
    
    cigarstring = 1S2M2I1M2D1M2S
    reference_start = 2
    >> r2q = reference2query(cigartuples, reference_start=2)
    (1, 2, 5, None, None, 6)
    """
    for cig_index in cigar_iterator(cigartuples, reference_start=reference_start):
        if cig_index.reference_index is not None:
            yield cig_index.query_index


def query2reference(cigartuples: CigarTuples, reference_start: int = 0) -> Iterator[Optional[int]]:
    """Create a generator the same size as the query
    that maps positions in the query to positions in the reference
    
    RPOS    0123  456789 # Index within the reference
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query 
    q2r =    N23NN4  7NN  
    
    cigarstring = 1S2M2I1M2D1M2S
    reference_start = 2
    >> q2r = query2reference(cigartuples, reference_start=2)
    [None, 2, 3, None, None, 4, 7, None, None]
    """
    for cig_index in cigar_iterator(cigartuples, reference_start=reference_start):
        if cig_index.query_index is not None:
            yield cig_index.reference_index


def query2cigar(cigartuples: CigarTuples, reference_start: int = 0) -> Iterator[Tuple[int, int]]:
    """Create a generator the same size as the query
    that maps positions in the query to CIGAR positions.
    
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query 
    
    CIGIND   011223  566 # Index of the cigar block
    CBLKIND  001010  001   # Index within the cigar block
    
    cigarstring = 1S2M2I1M2D1M2S
    reference_start = 2
    >> q2c = query2cigar(cigartuples, reference_start=2)
    [(0,0), (1,0), (1,1), (2,0), (2,1), (3,0), (5,0), (6,0), (6,1)]
    """
    for cig_index in cigar_iterator(cigartuples, reference_start=reference_start):
        if cig_index.query_index is not None:
            yield (cig_index.cigar_index, cig_index.cigar_block_index)

            
    
# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
