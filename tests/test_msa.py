"""Test multiple sequence operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import cigarstr2tup

import cigarmath as cm
from cigarmath import defn


def test_msa2cigartuples_basic():
    
    ref   = 'AAAAGACCCCCGACTAGCTAGCATGCT----ATCTAGCTAGCA'
    query = '----AACCCCCGAC----TAGCATGCTTTTTATCTAGCT----'
    # CIG        MMMMMMMMMMDDDDMMMMMMMMIIIIIMMMMMMMM
    
    cor_ref_start = 4
    cor_cigartuples = [(0, 10), (2, 4), (0, 9), (1, 4), (0, 8)]
    
    ref_start, cigartuples = cm.msa2cigartuples(ref, query)
    
    assert ref_start == cor_ref_start
    assert cigartuples == cor_cigartuples
    
    
def test_msa2cigartuples_extra_items():
    
    ref   = '--AAAAGACCCCCGACTA--GCTAGCATGCT----ATCTAGCTAGCA'
    query = '------AACCCCCGAC------TAGCATGCTTTTTATCTAGCT----'
    # CIG        MMMMMMMMMMDDDDMMMMMMMMIIIIIMMMMMMMM
    
    cor_ref_start = 4
    cor_cigartuples = [(0, 10), (2, 4), (0, 9), (1, 4), (0, 8)]
    
    ref_start, cigartuples = cm.msa2cigartuples(ref, query)
    
    assert ref_start == cor_ref_start
    assert cigartuples == cor_cigartuples
    
    
    
def test_msa2cigartuples_pre_items():
    
    ref   = '----AAAAGACCCCCGACTA--GCTAGCATGCT----ATCTAGCTAGCA---'
    query = 'TT------AACCCCCGAC------TAGCATGCTTTTTATCTAGCT----TTT'
    # CIG        MMMMMMMMMMDDDDMMMMMMMMIIIIIMMMMMMMM
    
    cor_ref_start = 4
    cor_cigartuples = [(4, 2), (0, 10), (2, 4), (0, 9), (1, 4), (0, 8), (4, 3)]
    
    ref_start, cigartuples = cm.msa2cigartuples(ref, query)
    
    assert ref_start == cor_ref_start
    assert cigartuples == cor_cigartuples
    
    
def test_softclipify():
    
    # REF    --AAAAGACCCCCGACTCGTTA---
    # QUE    TT----AACCCCCGAC----TAGCA
    # CIG    IIDDDDMMMMMMMMMMDDDDMMIII 
    # OUT    SS    MMMMMMMMMMDDDDMMSSS required_mapping = 1
    
    cigartuples = [(1, 2), (2,4), (0, 10), (2, 4), (0, 2), (1, 3)]
    
    cortuples = [(4, 2), (0, 10), (2, 4), (0, 2), (4, 3)]
    coroffset = 4
    
    newtuples, newoffset = cm.softclipify(cigartuples)
    
    assert newtuples == cortuples
    assert newoffset == coroffset
    
    
def test_decide_softclip_end():
    
    cigartuples = [(1, 2), (2,4), (0, 10), (2, 4), (0, 2), (1, 3)]
    
    left_cor_soft_size = 2
    left_cor_ind = 2
    
    right_cor_soft_size = 3
    right_cor_ind = 1
    
    left_ind, left_soft_sz = cm.clipping._decide_softclip_end(cigartuples, 1)
    
    assert left_cor_soft_size == left_soft_sz
    assert left_cor_ind == left_ind
    
    right_ind, right_soft_sz = cm.clipping._decide_softclip_end(cigartuples[::-1], 1)

    assert right_cor_soft_size == right_soft_sz
    assert right_cor_ind == right_ind
    
    slc = slice(left_ind, -right_ind or None)
    
    assert cigartuples[slc] == [(0, 10), (2, 4), (0, 2)]