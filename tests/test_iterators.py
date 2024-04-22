__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"


import cigarmath as cm
from cigarmath.defn import cigarstr2tup
from cigarmath.iterators import CigarIndex


def make_example():
    
    cigar = '1S2M2I1M2D1M2S'
    reference_start = 2
    
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
    
    return cm.cigarstr2tup(cigar), reference_start, correct

def check_index_list(guess_items, correct_items):

    for guess, cor in zip(guess_items, correct_items):
        print('cor  ', cor)
        print('guess', guess)
        print('\n\n')
        assert guess == cor
        
    assert len(guess_items) == len(correct_items)
    


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
    
    cigartups, ref_start, correct_indexes = make_example()
    cigiters = list(cm.cigar_iterator(cigartups, reference_start=ref_start))
    
    check_index_list(cigiters, correct_indexes)
        
        
def test_cigar_iterator_reference_slice():
    
    cigartups, ref_start, correct_indexes = make_example()
    
    cigiters = list(cm.cigar_iterator_reference_slice(cigartups,
                                                              reference_start=ref_start))
    print('Checking no boundary')
    check_index_list(cigiters, correct_indexes[1:])
    
    cigiters = list(cm.cigar_iterator_reference_slice(cigartups,
                                                              reference_start=ref_start,
                                                              region_reference_start = 5))
    print('Checking lower bound')
    check_index_list(cigiters, correct_indexes[6:])
    
    cigiters = list(cm.cigar_iterator_reference_slice(cigartups,
                                                              reference_start=ref_start,
                                                              region_reference_end = 5))
    print('Checking upper bound')
    check_index_list(cigiters, correct_indexes[1:6])
    
    
    cigiters = list(cm.cigar_iterator_reference_slice(cigartups,
                                                              reference_start=ref_start,
                                                              region_reference_start = 4,
                                                              region_reference_end = 7))
    print('Checking both bound')
    check_index_list(cigiters, correct_indexes[5:-3])
    
    
def test_liftover():
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
    
    cigartups, ref_start, _ = make_example()
    
    lift_ref_pos = [3, 4, 5, 6, 7]
    
    # right
    lift_que_pos = [2, 5, 6, 6, 6]
    lift_pos = list(cm.liftover(cigartups, *lift_ref_pos, reference_start=ref_start, edge='right'))
    
    assert lift_pos == lift_que_pos
    
    
    # left
    lift_que_pos = [2, 5, 5, 5, 6]
    lift_pos = list(cm.liftover(cigartups, *lift_ref_pos, reference_start=ref_start, edge='left'))
    
    assert lift_pos == lift_que_pos
    
    # None
    lift_que_pos = [2, 5, None, None, 6]
    lift_pos = list(cm.liftover(cigartups, *lift_ref_pos, reference_start=ref_start, edge=None))
    
    assert lift_pos == lift_que_pos