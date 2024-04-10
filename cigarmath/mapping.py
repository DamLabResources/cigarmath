"""Create mappings between read & reference based on CIGARs"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from dataclasses import dataclass

from cigarmath.defn import CONSUMES_REFERENCE, CONSUMES_QUERY, CLIPPING, BAM_CSOFT_CLIP


def cigar_iterator(cigartuples, reference_start=0):
    """Yields alignment, query, reference, and cigar indexes in alignment order.
    
    ALNPOS   01234567890 # Index of the entire alignment
    RPOS    0123  456789 # Index within the reference
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    CIGIND   01122344566 # Index of the cigar block
    CBLKIND  001010010   # Index within the cigar block
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query  
    
    
    for cigar_index in cigar_iterator(cigartuples, reference_start = 2):
        print(cigar_index)
    """
    
    if cigartuples[0][0] in CLIPPING:
        op, sz = cigartuples[0]
        yield from _left_clip_iterator(sz, cigar_op = op)
        
        alignment_index = sz-1
        query_index = sz-1
        cigartuples = cigartuples[1:]
        cigar_op_start=1
    else:
        query_index = -1
        alignment_index = -1
        cigar_op_start=0
        
    reference_index = reference_start-1
    
    for cigar_index, (cigar_op, size) in enumerate(cigartuples, cigar_op_start):
        
        # Precalc the changes for the block
        reference_delta = int(cigar_op in CONSUMES_REFERENCE)
        query_delta = int(cigar_op in CONSUMES_QUERY)
        
        for cigar_block_index in range(size):
            
            alignment_index += 1
            query_index += query_delta
            reference_index += reference_delta
            
            yield CigarIndex(alignment_index = alignment_index,
                             reference_index = reference_index if cigar_op in CONSUMES_REFERENCE else None,
                             query_index = query_index if cigar_op in CONSUMES_QUERY else None,
                             cigar_index = cigar_index,
                             cigar_block_index = cigar_block_index,
                             cigar_op = cigar_op)
            
            
def _left_clip_iterator(size, cigar_op = BAM_CSOFT_CLIP):
    "Handle left-clipping special case"
    
    for index in range(size):
        yield CigarIndex(alignment_index = index,
                         reference_index = None,
                         query_index = index,
                         cigar_index = 0,
                         cigar_block_index = index,
                         cigar_op = cigar_op)

@dataclass
class CigarIndex:
    alignment_index: int
    reference_index: int
    query_index: int
    cigar_index: int
    cigar_block_index: int
    cigar_op: int
        
        

def reference2query(cigartuples, reference_start=0):
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
    
    
def query2reference(cigartuples, reference_start=0):
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
    
    
def query2cigar(cigartuples, reference_start=0):
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
