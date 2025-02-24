"""Unit tests for pileup functions"""

from cigarmath.pileup import depth
from cigarmath.defn import cigarstr2tup
from collections import Counter

def test_basic_depth():
    """Test basic pileup depth calculation"""
    # Example from docstring
    
    cigarstring = '4M 2I 2M 2D 3M'    
    #                 MMMMIIMMDDMMMM
    query_sequence = "AAAAACCGGCC"
    
    cigartuples = cigarstr2tup(cigarstring)
    
    result = depth(cigartuples, query_sequence=query_sequence)
    
    assert result[0] == Counter({'A': 1})
    assert result[1] == Counter({'A': 1})
    assert result[2] == Counter({'A': 1})
    assert result[3] == Counter({'A': 1})

    assert result[4] == Counter({'C': 1})  
    assert result[5] == Counter({'G': 1})
    
    # Deletion
    # Deletion
    
    # GCC
    assert result[8] == Counter({'G': 1})
    assert result[9] == Counter({'C': 1})
    assert result[10] == Counter({'C': 1})

def test_reference_start():
    """Test pileup with non-zero reference start"""
    cigartuples = [(0,3), (2,2), (0,2)]  # 3M2D2M
    query_sequence = "AAATT"
    
    result = depth(cigartuples, reference_start=10, query_sequence=query_sequence)
    
    assert result[10] == Counter({'A': 1})
    assert result[11] == Counter({'A': 1})
    assert result[12] == Counter({'A': 1})
    assert result[13] == Counter({'-': 1})
    assert result[14] == Counter({'-': 1})
    assert result[15] == Counter({'T': 1})
    assert result[16] == Counter({'T': 1})

def test_no_query_sequence():
    """Test pileup when query sequence is not provided"""
    cigartuples = [(0,3), (2,1), (0,2)]  # 3M1D2M
    
    result = depth(cigartuples)
    
    assert result[0] == Counter({'.': 1})
    assert result[1] == Counter({'.': 1})
    assert result[2] == Counter({'.': 1})
    assert result[3] == Counter({'-': 1})
    assert result[4] == Counter({'.': 1})
    assert result[5] == Counter({'.': 1})

def test_previous_count():
    """Test updating previous pileup counts"""
    # First alignment
    cigartuples1 = [(0,3)]  # 3M
    query_sequence1 = "AAA"
    result1 = depth(cigartuples1, query_sequence=query_sequence1)
    
    # Second alignment
    cigartuples2 = [(0,3)]  # 3M
    query_sequence2 = "TTT"
    result2 = depth(cigartuples2, query_sequence=query_sequence2, previous_count=result1)
    
    # Check combined counts
    assert result2[0] == Counter({'A': 1, 'T': 1})
    assert result2[1] == Counter({'A': 1, 'T': 1})
    assert result2[2] == Counter({'A': 1, 'T': 1})


def test_skipped_regions():
    """Test pileup with skipped regions (N in CIGAR)"""
    cigartuples = [(0,2), (3,3), (0,2)]  # 2M3N2M
    query_sequence = "AATT"
    
    result = depth(cigartuples, query_sequence=query_sequence)
    
    assert result[0] == Counter({'A': 1})
    assert result[1] == Counter({'A': 1})
    assert result[2] == Counter({'-': 1})
    assert result[3] == Counter({'-': 1})
    assert result[4] == Counter({'-': 1})
    assert result[5] == Counter({'T': 1})
    assert result[6] == Counter({'T': 1}) 