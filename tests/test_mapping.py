import cigarmath as cm
        
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
    
