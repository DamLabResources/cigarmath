"""Functions for combining CIGAR strings"""

from typing import Tuple, List, Optional
from functools import partial
from cigarmath.block import reference_block, query_block
from cigarmath.defn import CigarTuples, BAM_CDEL, CONSUMES_QUERY, CONSUMES_REFERENCE
from cigarmath.cigarmath import collapse_adjacent_blocks
from cigarmath.clipping import declip

def combine_adjacent_alignments(
    first: Tuple[int, CigarTuples],
    second: Tuple[int, CigarTuples]
) -> Tuple[int, CigarTuples]:
    """Combine two adjacent alignments, adding deletion if there's a gap or trimming if there's overlap.
    
    Args:
        first: Tuple of (reference_start, cigartuples) for first alignment
        second: Tuple of (reference_start, cigartuples) for second alignment
    
    Returns:
        Tuple of (reference_start, combined_cigartuples)
        
    Example:
        First alignment: REF: AAAAGATC--CCC  (ref_start=4)
                         QRY: AAAA-ATCGGCCC
                         CGT: 4M1D3M2I3M
                         
        Second alignment:REF:        CCCTAG  (ref_start=12) 
                         QRY:         CCTAG
                         CGT:         3M2M
                         
        Combined:        REF: AAAAGATC--CCCCCTAG
                        QRY: AAAA-ATCGGCCCCCTAG  
                        CGT: 4M1D3M2I6M2M
    """
    first_start, first_cigars = first[0], declip(first[1])
    second_start, second_cigars = second[0], declip(second[1])
    
    assert first_start < second_start, 'Can only combine adjacent alignments'

    # Get the reference blocks to check for overlap/gap
    first_ref_start, first_ref_end = reference_block(first_cigars, first_start)
    
    # Calculate gap/overlap size
    gap_size = second_start - first_ref_end
    
    if gap_size > 0:
        # There's a gap - add deletion operation between alignments
        combined_cigars = list(first_cigars)
        combined_cigars.append((BAM_CDEL, gap_size))
        combined_cigars.extend(second_cigars)
        return first_start, tuple(collapse_adjacent_blocks(combined_cigars))
        
    elif gap_size < 0:
        # There's an overlap - need to trim second alignment
        # Get query coordinates to determine where to trim
        first_q_start, first_q_end = query_block(first_cigars)
        second_q_start, second_q_end = query_block(second_cigars)
        
        # Calculate how many query bases to trim from second alignment
        # This is based on how many reference bases overlap
        overlap = -gap_size
        
        # Trim the overlapping portion from the start of second alignment
        trimmed_start, trimmed_cigars = trim_alignment(
            second_start,
            second_cigars,
            left=overlap
        )
        print('second cigars', second_cigars, 'trimming', overlap)
        print('trimmed cigars', trimmed_cigars)
        
        # Combine the alignments
        combined_cigars = list(first_cigars)
        combined_cigars.extend(trimmed_cigars)
        return first_start, tuple(collapse_adjacent_blocks(combined_cigars))
        
    else:
        # Alignments are perfectly adjacent
        combined_cigars = list(first_cigars)
        combined_cigars.extend(second_cigars)
        return first_start, tuple(collapse_adjacent_blocks(combined_cigars))

def _trim(
    cigartuples: CigarTuples,
    trim_amount: int,
    from_start: bool = True
) -> Tuple[CigarTuples, int]:
    """Helper function to trim query bases from either end of alignment.
    
    Args:
        cigartuples: List of (operation, length) tuples
        trim_amount: Number of query bases to trim
        from_start: If True, trim from start, if False trim from end
        
    Returns:
        Tuple of (trimmed_cigartuples, reference_position_delta)
    """
    if trim_amount == 0:
        return cigartuples, 0
        
    new_tuples = []
    ref_pos_delta = 0
    to_trim = trim_amount
    
    # Reverse tuples if trimming from end
    tuples = reversed(cigartuples) if not from_start else cigartuples
    
    # Get the appropriate method for appending or inserting
    extend_func = partial(new_tuples.append) if from_start else partial(new_tuples.insert, 0)
    
    for op, length in tuples:
        if to_trim > 0:
            if op in CONSUMES_QUERY:
                if length <= to_trim:
                    # Skip this operation entirely
                    to_trim -= length
                    if op in CONSUMES_REFERENCE:
                        ref_pos_delta += length
                    continue
                else:
                    # Partial operation
                    new_length = length - to_trim
                    extend_func((op, new_length))
                    
                    if op in CONSUMES_REFERENCE:
                        ref_pos_delta += to_trim
                    to_trim = 0
            else:
                # Non-query consuming op
                if op in CONSUMES_REFERENCE:
                    ref_pos_delta += length
        else:
            # Past trimming
            extend_func((op, length))
                
    return tuple(new_tuples), ref_pos_delta


def trim_alignment(
    ref_start: int,
    cigartuples: CigarTuples,
    left: int = 0,
    right: int = 0,
    add_clipping: Optional[str] = None
) -> Tuple[int, CigarTuples]:
    """Trim query bases from left and/or right of alignment.
    
    Args:
        ref_start: Starting position on reference
        cigartuples: List of (operation, length) tuples
        left: Number of query bases to trim from left
        right: Number of query bases to trim from right
        add_clipping: If 'soft' or 'hard', add clipping operations to maintain query length
        
    Returns:
        Tuple of (new_ref_start, trimmed_cigartuples)
        
    Example:
        ref_start = 10
        cigartuples = [(0,5), (2,1), (0,3)]  # 5M1D3M
        left = 2
        right = 1
        add_clipping = 'soft'
        
        Returns: (12, [(4,2), (0,3), (2,1), (0,2), (4,1)])  # 2S3M1D2M1S
    """
    if (left == 0) and (right == 0):
        return ref_start, cigartuples
        
    if add_clipping not in (None, 'soft', 'hard'):
        raise ValueError("add_clipping must be None, 'soft', or 'hard'")
        
    # Handle left trimming first
    trimmed_cigars, ref_delta = _trim(cigartuples, left, from_start=True)
    
    # Then right trimming if needed
    if right > 0:
        trimmed_cigars, _ = _trim(trimmed_cigars, right, from_start=False)
        
    # Add clipping if requested
    if add_clipping:
        from cigarmath.defn import BAM_CSOFT_CLIP, BAM_CHARD_CLIP
        clip_op = BAM_CSOFT_CLIP if add_clipping == 'soft' else BAM_CHARD_CLIP
        
        result = []
        if left > 0:
            result.append((clip_op, left))
        result.extend(trimmed_cigars)
        if right > 0:
            result.append((clip_op, right))
        trimmed_cigars = tuple(result)
        
    return ref_start + ref_delta, trimmed_cigars

def combine_multiple_alignments(
    alignments: List[Tuple[int, CigarTuples]],
    allowed_overlap: int = 0
) -> Tuple[int, CigarTuples]:
    """Combine multiple alignments after sorting by query position and validating.
    
    Args:
        alignments: List of (reference_start, cigartuples) tuples
        allowed_overlap: Maximum number of query bases that can overlap between adjacent alignments
        
    Returns:
        Tuple of (reference_start, combined_cigartuples)
        
    Raises:
        ValueError: If alignments overlap more than allowed or are not sequential
        
    Example:
        alignments = [
            (10, [(0,5), (2,2), (0,3)]),  # 5M2D3M at ref:10
            (25, [(0,4), (1,2), (0,2)]),  # 4M2I2M at ref:25
            (35, [(0,5)])                  # 5M at ref:35
        ]
        allowed_overlap = 2
        
        Returns: (10, [(0,5), (2,2), (0,7), (1,2), (0,7)])
    """
    if not alignments:
        raise ValueError("No alignments provided")
    if len(alignments) == 1:
        return alignments[0]
        
    # Sort alignments by query start position
    sorted_alns = sorted(
        alignments,
        key=lambda x: query_block(x[1])[0]
    )
    
    # Validate that alignments don't overlap more than allowed
    # and that they are sequential in reference space
    prev_query_end = None
    prev_ref_end = None
    
    for ref_start, cigars in sorted_alns:
        q_start, q_end = query_block(cigars)
        r_start, r_end = reference_block(cigars, ref_start)
        
        if prev_query_end is not None:
            overlap = prev_query_end - q_start
            if (prev_query_end - q_start) > allowed_overlap:
                raise ValueError(
                    f"Alignments overlap by {overlap} bases in query space "
                    f"(maximum allowed: {allowed_overlap}): "
                    f"Previous alignment ends at {prev_query_end}, "
                    f"but next alignment starts at {q_start}"
                )
            if (prev_ref_end - r_start) > allowed_overlap:
                raise ValueError(
                    f"Alignments are not sequential in reference space: "
                    f"Previous alignment ends at {prev_ref_end}, "
                    f"but next alignment starts at {r_start}"
                )
                
        prev_query_end = q_end
        prev_ref_end = r_end
    
    # Combine alignments sequentially
    result = sorted_alns[0]
    for next_aln in sorted_alns[1:]:
        result = combine_adjacent_alignments(result, next_aln)
        
    return result
