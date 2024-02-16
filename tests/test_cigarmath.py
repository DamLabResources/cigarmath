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


def test_left_clipping():
    "Test detecting left clipping"

    cigar = "30M"
    guess = cigarmath.left_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected no clipping"

    cigar = "20S30M10S"
    guess = cigarmath.left_clipping(cigarstr2tup(cigar))
    assert guess == 20, "Expected 20 left-clipping"

    cigar = "30M10S"
    guess = cigarmath.left_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected only right-clipping"

    cigar = "20H30M10H"
    guess = cigarmath.left_clipping(cigarstr2tup(cigar))
    assert guess == 20, "Expected 20 left-hard-clipping"

    cigar = "20H30M10H"
    guess = cigarmath.left_clipping(cigarstr2tup(cigar), with_hard=False)
    assert guess == 0, "Expected excluding of left-hard-clipping"


def test_right_clipping():
    "Test detecting right clipping"

    cigar = "30M"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected no clipping"

    cigar = "20S30M10S"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 10 right-clipping"

    cigar = "30M10S"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 10 right-clipping"

    cigar = "10S30M"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar))
    assert guess == 0, "Expected only left-clipping"

    cigar = "20H30M10H"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar))
    assert guess == 10, "Expected 20 left-hard-clipping"

    cigar = "20H30M10H"
    guess = cigarmath.right_clipping(cigarstr2tup(cigar), with_hard=False)
    assert guess == 0, "Expected excluding of right-hard-clipping"


def test_reference_block():
    "Test determining the aligned reference block"

    ref_start = 10

    cigar = "30M"
    guess = cigarmath.reference_block(
        cigarstr2tup(cigar), reference_start=ref_start
    )
    assert (ref_start, ref_start + 30) == guess

    cigar = "20S30M10S"
    guess = cigarmath.reference_block(
        cigarstr2tup(cigar), reference_start=ref_start
    )
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Soft-Clipping should not alter the reference block"

    cigar = "20H30M10H"
    guess = cigarmath.reference_block(
        cigarstr2tup(cigar), reference_start=ref_start
    )
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Hard-Clipping should not alter the reference block"

    cigar = "20S25M10I5M10S"
    guess = cigarmath.reference_block(
        cigarstr2tup(cigar), reference_start=ref_start
    )
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Insertions should not alter the reference block"

    cigar = "20S25M10D5M10S"
    guess = cigarmath.reference_block(
        cigarstr2tup(cigar), reference_start=ref_start
    )
    assert (
        ref_start,
        ref_start + 40,
    ) == guess, "Deletions should alter the reference block"


def test_query_block():
    "Test determing the query alignment block"

    cigar = "30M"
    guess = cigarmath.query_block(cigarstr2tup(cigar))
    assert (0, 30) == guess

    cigar = "20S30M10S"
    guess = cigarmath.query_block(cigarstr2tup(cigar))
    assert (20, 50) == guess, "Soft-Clipping should shift the query block"

    cigar = "20H30M10H"
    guess = cigarmath.query_block(cigarstr2tup(cigar))
    assert (20, 50) == guess, "Hard-Clipping should shift the query block"

    cigar = "20S25M10I5M10S"
    guess = cigarmath.query_block(cigarstr2tup(cigar))
    assert (20, 20 + 25 + 10 + 5) == guess, "Insertions should shift the query block"

    cigar = "20S25M10D5M10S"
    guess = cigarmath.query_block(cigarstr2tup(cigar))
    assert (20, 20 + 25 + 5) == guess, "Deletions should not shift the query block"


def test_declip():
    "Test declipping cigartuples"

    cigar = "30M"
    cigartups = cigarstr2tup(cigar)
    guess = cigarmath.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30S30M20S"
    cigartups = cigarstr2tup(cigar)
    guess = cigarmath.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30S30M"
    cigartups = cigarstr2tup(cigar)
    guess = cigarmath.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)

    cigar = "30M20S"
    cigartups = cigarstr2tup(cigar)
    guess = cigarmath.declip(cigartups)
    correct = [(defn.BAM_CMATCH, 30)]
    check_cigartuples(guess, correct)


def test_block_overlap():
    "Test detecting overlapping blocks"

    tests = [
        ((10, 15), (12, 16), 3),  # One edge overlaps
        ((10, 15), (15, 20), 0),  # adjacent
        ((10, 20), (15, 17), 2),  # fully encased
        ((10, 15), (17, 20), -2),  # non-overlapping
    ]

    for block_a, block_b, correct in tests:
        guess = cigarmath.block_overlap_length(block_a, block_b)
        assert (
            correct == guess
        ), f"Expected {correct} overlap but got {guess} for {block_a}, {block_b}"

        guess = cigarmath.block_overlap_length(block_b, block_a)
        assert (
            correct == guess
        ), f"Expected {correct} overlap but got {guess} for {block_b}, {block_a}"


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


def test_is_hard_clipping():
    "Testing detecting hard-clipping"

    hard = map(cigarstr2tup, ["3H2M", "3M2H"])
    nothard = map(cigarstr2tup, ["2M", "3M2S"])

    assert all(cigarmath.is_hard_clipped(h) for h in hard)
    assert not any(cigarmath.is_hard_clipped(n) for n in nothard)
