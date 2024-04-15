"""A set of iterators over CIGARs"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from dataclasses import dataclass

from cigarmath.defn import CONSUMES_REFERENCE, CONSUMES_QUERY, CLIPPING, BAM_CSOFT_CLIP


@dataclass
class CigarIndex:
    "A helper class for holding information about each position in the alignment"
    
    alignment_index: int
    reference_index: int
    query_index: int
    cigar_index: int
    cigar_block_index: int
    cigar_op: int


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


        

def cigar_iterator_reference_slice(cigartuples, reference_start=0, region_reference_start=None, region_reference_end=None):
    "Return only a slice of the cigar_iterator based on a reference region"
    
    started = False
    
    for cigar_index in cigar_iterator(cigartuples, reference_start=reference_start):
        
        if (cigar_index.reference_index is None) and (started == False):
            # Clipped or insertion BEFORE starting
            # so skip
            pass
        elif (cigar_index.reference_index is None) and started:
            # Insertion within the region, so yielding
            yield cigar_index
        elif (region_reference_start is not None) and (cigar_index.reference_index < region_reference_start):
            # Has a lower-bound, this is below it 
            pass
        elif (region_reference_end is not None) and (cigar_index.reference_index >= region_reference_end):
            # Has an upper-bound, this is above it 
            break
        else:
            started = True
            yield cigar_index




# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
