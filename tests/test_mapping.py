__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"


import cigarmath as cm
from cigarmath.defn import cigarstr2tup
from cigarmath.mapping import CigarIndex


def test_cigar_iterator():
    """
    
    ALNPOS   01234567890 # Index of the entire alignment
    RPOS    0123  456789 # Index within the reference
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    CIGIND   01122344566 # Index of the cigar block
    CBLKIND  001010010   # Index within the cigar block
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query 
    """
    
    cigar = '1S2M2I1M2D1M2S'
    
    
    correct = [CigarIndex(alignment_index = 0,
                          reference_index = None,
                          query_index = 0,
                          cigar_index = 0,
                          cigar_block_index=0,
                          cigar_op = 4),
               CigarIndex(alignment_index = 1,
                          reference_index = 2,
                          query_index = 1,
                          cigar_index = 1,
                          cigar_block_index=0,
                          cigar_op = 0),
               CigarIndex(alignment_index = 2,
                          reference_index = 3,
                          query_index = 2,
                          cigar_index = 1,
                          cigar_block_index=1,
                          cigar_op = 0),
               CigarIndex(alignment_index = 3,
                          reference_index = None,
                          query_index = 3,
                          cigar_index = 2,
                          cigar_block_index=0,
                          cigar_op = 1),
               CigarIndex(alignment_index = 4,
                          reference_index = None,
                          query_index = 4,
                          cigar_index = 2,
                          cigar_block_index=1,
                          cigar_op = 1),
               CigarIndex(alignment_index = 5,
                          reference_index = 4,
                          query_index = 5,
                          cigar_index = 3,
                          cigar_block_index=0,
                          cigar_op = 0),
               CigarIndex(alignment_index = 6,
                          reference_index = 5,
                          query_index = None,
                          cigar_index = 4,
                          cigar_block_index=0,
                          cigar_op = 2),
               CigarIndex(alignment_index = 7,
                          reference_index = 6,
                          query_index = None,
                          cigar_index = 4,
                          cigar_block_index=1,
                          cigar_op = 2),
               CigarIndex(alignment_index = 8,
                          reference_index = 7,
                          query_index = 6,
                          cigar_index = 5,
                          cigar_block_index=0,
                          cigar_op = 0),
               
               CigarIndex(alignment_index = 9,
                          reference_index = None,
                          query_index = 7,
                          cigar_index = 6,
                          cigar_block_index=0,
                          cigar_op = 4),
               
               CigarIndex(alignment_index = 10,
                          reference_index = None,
                          query_index = 8,
                          cigar_index = 6,
                          cigar_block_index=1,
                          cigar_op = 4),
               ]
    
    
    cigartups = cm.cigarstr2tup(cigar)
    
    cigiters = list(cm.mapping.cigar_iterator(cigartups, reference_start=2))
    
    assert len(correct) == len(cigiters)
    
    for guess, cor in zip(cigiters, correct):
        print('cor  ', cor)
        print('guess', guess)
        print('\n\n')
        assert guess == cor
        
        
def test_reference2query():
    """RPOS    0123  456789 # Index within the reference
    REF     AAGA--CTTCGG
    CIGAR    SMMIIMDDMSS
    QRY     -xAAGGC--Cxx
    QPOS     012345  678 # Index within the query 
    
    r2q =     12  5NN6  
    
    cigarstring = 1S2M2I1M2D1M2S
    reference_start = 2
    >> r2q = reference2query(cigartuples, reference_start=2)
    [1, 2, 5, None, None, 6]
    """
    
    cigar = '1S2M2I1M2D1M2S'
    cigartuples = cm.cigarstr2tup(cigar)
    r2q = list(cm.reference2query(cigartuples, reference_start=2))
    correct = [1, 2, 5, None, None, 6]
    
    assert r2q == correct
    
    
def test_query2reference():
    """
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
    
    cigar = '1S2M2I1M2D1M2S'
    cigartuples = cm.cigarstr2tup(cigar)
    r2q = list(cm.query2reference(cigartuples, reference_start=2))
    correct = [None, 2, 3, None, None, 4, 7, None, None]
    
    assert r2q == correct