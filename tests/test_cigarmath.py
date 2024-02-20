"""Test basic cigarmath operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath import cigarmath
from cigarmath import defn

from cigarmath.defn import cigarstr2tup


def check_cigartuples(guess, correct):
    "Compare two lists of cigartuples"
    assert len(guess) == len(correct)
    for (gop, gsz), (cop, csz) in zip(guess, correct):
        assert gop == cop
        assert gsz == csz


def test_cigarstr2tuples():
    "Test converting cigarstrings to cigartuples"

    cigar = "30M"
    guess = cigarstr2tup(cigar)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30S30M20S"
    guess = cigarstr2tup(cigar)
    correct = [
        (defn.BAM_CSOFT_CLIP, 30),
        (defn.BAM_CMATCH, 30),
        (defn.BAM_CSOFT_CLIP, 20),
    ]
    check_cigartuples(guess, correct)





def test_simplify_blocks():
    "Test collapsing adjacent cigar blocks of the same type"

    # Test safe passthrough
    cigar = "30M20I"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.simplify_blocks(cigartups))
    correct = [
        (defn.BAM_CMATCH, 30),
        (defn.BAM_CINS, 20),
    ]
    check_cigartuples(guess, correct)

    # Test replace without collapse
    cigar = "15=20X15I"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.simplify_blocks(cigartups, collapse=False))
    correct = [
        (defn.BAM_CMATCH, 15),
        (defn.BAM_CMATCH, 20),
        (defn.BAM_CINS, 15),
    ]
    check_cigartuples(guess, correct)

    # Test replace with collapse
    guess = list(cigarmath.simplify_blocks(cigartups, collapse=True))
    correct = [
        (defn.BAM_CMATCH, 35),
        (defn.BAM_CINS, 15),
    ]
    check_cigartuples(guess, correct)


def test_collapse_adjacent_blocks():
    "Test collapsing adjacent cigar blocks of the same type"

    # Test safe passthrough
    cigar = "30M"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.collapse_adjacent_blocks(cigartups))
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    # Test replace without collapse
    cigar = "15M20M15I"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.collapse_adjacent_blocks(cigartups))
    correct = [
        (defn.BAM_CMATCH, 35),
        (defn.BAM_CINS, 15),
    ]
    check_cigartuples(guess, correct)
