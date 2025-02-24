"""Functions for calculating pileup statistics from CIGAR strings"""

from collections import defaultdict, Counter
from typing import Dict, Optional
from cigarmath.defn import (
    CigarTuples,
    CONSUMES_REFERENCE,
    CONSUMES_QUERY,
    BAM_CDEL,
    BAM_CREF_SKIP,
)

def depth(
    cigartuples: CigarTuples,
    reference_start: int = 0,
    query_sequence: Optional[str] = None,
    previous_count: Optional[Dict[int, Counter]] = None,
) -> Dict[int, Counter]:
    """Calculate the pileup depth at each reference position.
    
    Args:
        cigartuples: List of CIGAR operations and their lengths
        reference_start: Starting position on reference (default: 0)
        query_sequence: Query sequence string (default: None, will use '.' for bases)
        previous_count: Previous pileup counts to update (default: None)
    
    Returns:
        Dict mapping reference positions to Counter of bases at that position
        
    Example:
        REF:     AAAAGACC--CCC
        QRY:     AAAA-ACCGGCCC
        CIGAR: 4M1D3M2I3M
        
        >>> depth([(0,4), (2,1), (0,3), (1,2), (0,3)], query_sequence="AAAAACCGGCCC")
        {
            0: Counter({'A': 1}),
            1: Counter({'A': 1}),
            2: Counter({'A': 1}),
            3: Counter({'A': 1}),
            4: Counter({'-': 1}),  # Deletion
            5: Counter({'A': 1}),
            6: Counter({'C': 1}),
            7: Counter({'C': 1}),
            8: Counter({'C': 1}),
            9: Counter({'C': 1}),
        }
    """
    # Initialize counts or use previous
    counts = defaultdict(Counter)
    if previous_count is not None:
        counts.update(previous_count)
        
    # Track current positions
    ref_pos = reference_start
    query_pos = 0
    
    for op, length in cigartuples:
        # Handle operations that consume reference
        if op in CONSUMES_REFERENCE:
            for i in range(length):
                if op == BAM_CDEL or op == BAM_CREF_SKIP:
                    # Mark deletions with '-'
                    counts[ref_pos]['-'] += 1
                else:
                    # Get base from query sequence or use '.' if not provided
                    base = '.'
                    if query_sequence is not None and query_pos < len(query_sequence):
                        base = query_sequence[query_pos]
                    counts[ref_pos][base] += 1
                    if op in CONSUMES_QUERY:
                        query_pos += 1
                ref_pos += 1
        # Handle insertions (don't advance ref_pos)
        elif op in CONSUMES_QUERY:
            query_pos += length
            
    return counts 