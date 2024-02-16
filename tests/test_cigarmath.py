"""Test basic cigarmath operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath import cigarmath
from cigarmath import defn

from cigarmath.defn import cigarstring_to_cigartuples as cigarstr2tup

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












def test_reference_deletion_blocks():
    "Test detecting large deletions in cigartuples"

    cigar = "30M10D20S"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.reference_deletion_blocks(cigartups))
    correct = [(30, 40)]
    assert guess == correct

    guess = list(cigarmath.reference_deletion_blocks(cigartups, reference_start=10))
    correct = [(40, 50)]
    assert guess == correct

    guess = list(cigarmath.reference_deletion_blocks(cigartups, min_size=20))
    correct = []
    assert guess == correct

    cigar = "30M10N30M100D10M10I10M50N10M"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.reference_deletion_blocks(cigartups))
    correct = [(30, 40), (70, 170), (190, 240)]
    assert guess == correct

    guess = list(cigarmath.reference_deletion_blocks(cigartups, min_size=20))
    correct = [(70, 170), (190, 240)]
    assert guess == correct

    # Test handles no-deletions
    cigar = "300M"
    cigartups = cigarstr2tup(cigar)
    guess = list(cigarmath.reference_deletion_blocks(cigartups))
    correct = []
    assert guess == correct


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


