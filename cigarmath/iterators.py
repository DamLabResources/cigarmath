"""A set of iterators over CIGARs"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from dataclasses import dataclass
from itertools import dropwhile
from functools import partial
from typing import Optional

from cigarmath.defn import CONSUMES_REFERENCE, CONSUMES_QUERY, CLIPPING, BAM_CSOFT_CLIP
from cigarmath.block import reference_offset


@dataclass
class CigarIndex:
    "A helper class for holding information about each position in the alignment"
    
    alignment_index: int
    reference_index: int
    query_index: int
    cigar_index: int
    cigar_block_index: int
    cigar_op: int
    
    query_letter: Optional[str] = None
    query_quality: Optional[int] = None
    reference_letter: Optional[str] = None


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


def iterator_attach(cigar_index_iterator, reference_sequence=None, query_sequence=None, query_qualities=None):
    "Attach reference and query information to a cigar_index_iterator"
    
    for cigar_index in cigar_index_iterator:
        if reference_sequence and (cigar_index.reference_index is not None):
            cigar_index.reference_letter = reference_sequence[cigar_index.reference_index]
        if query_sequence and (cigar_index.query_index is not None):
            cigar_index.query_letter = query_sequence[cigar_index.query_index]
            if query_qualities:
                cigar_index.query_quality = query_qualities[cigar_index.query_index]
        yield cigar_index
            
        
    
            
            
def is_sorted_ascending(lst):
    return all(lst[i] <= lst[i+1] for i in range(len(lst) - 1))
            
def _before_site_on_ref(site, index):
    "return True if this index if before this site on the reference"
    
    if index.reference_index is None:
        return True
    return index.reference_index < site-1
    
            
def liftover(cigartuples, *reference_sites, reference_start=0, edge = 'left', small_indel_limit = 10):
    "Create an iterator which yields the query index for each reference site"
    
    assert edge in {'left', 'right', None}
    
    cigar_indexes = list(cigar_iterator(cigartuples, reference_start=reference_start))
    
    for site in reference_sites:
        yield _liftover_index(cigar_indexes, site, edge, small_indel_limit)
    
    
def _liftover_index(cigar_indexes, site, edge, limit):
    
    for n in range(len(cigar_indexes)):
        
        if cigar_indexes[n].reference_index == site:
            if (cigar_indexes[n].query_index is not None) or (edge is None):
                # Mapped or don't care
                return cigar_indexes[n].query_index
            elif (cigar_indexes[n].query_index is None) and (edge == 'left'):
                # Check backwards along the alignment to find a query position that maps
                for left_check in range(n-1, max(n-limit, 0), -1):
                    if cigar_indexes[left_check].query_index is not None:
                        return cigar_indexes[left_check].query_index
                # Checked back past the limit, must be a long deletion
                return None
            elif (cigar_indexes[n].query_index is None) and (edge == 'right'):
                # Check forwards along the alignment to find a query position that maps
                for right_check in range(n+1, min(n+limit+1, len(cigar_indexes))):
                    if cigar_indexes[right_check].query_index is not None:
                        return cigar_indexes[right_check].query_index
                # Checked back past the limit, must be a long deletion
                return None


# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
