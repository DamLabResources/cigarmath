"""Test basic clipping operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import cigarstr2tup

import cigarmath as cm
from cigarmath import defn


def check_cigartuples(guess, correct):
    "Compare two lists of cigartuples"
    assert len(guess) == len(correct)
    for (gop, gsz), (cop, csz) in zip(guess, correct):
        assert gop == cop
        assert gsz == csz


def test_left_clipping():
    "Test detecting left clipping"

    cigar = "30M"
    guess = cm.left_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected no clipping"

    cigar = "20S30M10S"
    guess = cm.left_clipping(cigarstr2tup(cigar))
    assert guess == 20, "Expected 20 left-clipping"

    cigar = "30M10S"
    guess = cm.left_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected only right-clipping"

    cigar = "20H30M10H"
    guess = cm.left_clipping(cigarstr2tup(cigar))
    assert guess == 20, "Expected 20 left-hard-clipping"

    cigar = "20H30M10H"
    guess = cm.left_clipping(cigarstr2tup(cigar), with_hard=False)
    assert guess == 0, "Expected excluding of left-hard-clipping"


def test_right_clipping():
    "Test detecting right clipping"

    cigar = "30M"
    guess = cm.right_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected no clipping"

    cigar = "20S30M10S"
    guess = cm.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 10 right-clipping"

    cigar = "30M10S"
    guess = cm.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 10 right-clipping"

    cigar = "10S30M"
    guess = cm.right_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected only left-clipping"

    cigar = "20H30M10H"
    guess = cm.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 20 left-hard-clipping"

    cigar = "20H30M10H"
    guess = cm.right_clipping(cigarstr2tup(cigar), with_hard=False)
    assert guess == 0, "Expected excluding of right-hard-clipping"


def test_is_hard_clipping():
    "Testing detecting hard-clipping"

    hard = map(cigarstr2tup, ["3H2M", "3M2H"])
    nothard = map(cigarstr2tup, ["2M", "3M2S"])

    assert all(cm.is_hard_clipped(h) for h in hard)
    assert not any(cm.is_hard_clipped(n) for n in nothard)


def test_declip():
    "Test declipping cigartuples"

    cigar = "30M"
    cigartups = cigarstr2tup(cigar)
    guess = cm.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30S30M20S"
    cigartups = cigarstr2tup(cigar)
    guess = cm.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30S30M"
    cigartups = cigarstr2tup(cigar)
    guess = cm.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30M20S"
    cigartups = cigarstr2tup(cigar)
    guess = cm.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)
