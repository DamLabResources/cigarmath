import cigarmath as cm
from cigarmath.defn import cigarstring_to_cigartuples as cigarstr2tup


def test_reference_block():
    "Test determining the aligned reference block"

    ref_start = 10

    cigar = "30M"
    guess = cm.reference_block(cigarstr2tup(cigar), reference_start=ref_start)
    assert (ref_start, ref_start + 30) == guess

    cigar = "20S30M10S"
    guess = cm.reference_block(cigarstr2tup(cigar), reference_start=ref_start)
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Soft-Clipping should not alter the reference block"

    cigar = "20H30M10H"
    guess = cm.reference_block(cigarstr2tup(cigar), reference_start=ref_start)
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Hard-Clipping should not alter the reference block"

    cigar = "20S25M10I5M10S"
    guess = cm.reference_block(cigarstr2tup(cigar), reference_start=ref_start)
    assert (
        ref_start,
        ref_start + 30,
    ) == guess, "Insertions should not alter the reference block"

    cigar = "20S25M10D5M10S"
    guess = cm.reference_block(cigarstr2tup(cigar), reference_start=ref_start)
    assert (
        ref_start,
        ref_start + 40,
    ) == guess, "Deletions should alter the reference block"


def test_query_block():
    "Test determing the query alignment block"

    cigar = "30M"
    guess = cm.query_block(cigarstr2tup(cigar))
    assert (0, 30) == guess

    cigar = "20S30M10S"
    guess = cm.query_block(cigarstr2tup(cigar))
    assert (20, 50) == guess, "Soft-Clipping should shift the query block"

    cigar = "20H30M10H"
    guess = cm.query_block(cigarstr2tup(cigar))
    assert (20, 50) == guess, "Hard-Clipping should shift the query block"

    cigar = "20S25M10I5M10S"
    guess = cm.query_block(cigarstr2tup(cigar))
    assert (20, 20 + 25 + 10 + 5) == guess, "Insertions should shift the query block"

    cigar = "20S25M10D5M10S"
    guess = cm.query_block(cigarstr2tup(cigar))
    assert (20, 20 + 25 + 5) == guess, "Deletions should not shift the query block"


def test_block_overlap():
    "Test detecting overlapping blocks"

    tests = [
        ((10, 15), (12, 16), 3),  # One edge overlaps
        ((10, 15), (15, 20), 0),  # adjacent
        ((10, 20), (15, 17), 2),  # fully encased
        ((10, 15), (17, 20), -2),  # non-overlapping
    ]

    for block_a, block_b, correct in tests:
        guess = cm.block_overlap_length(block_a, block_b)
        assert (
            correct == guess
        ), f"Expected {correct} overlap but got {guess} for {block_a}, {block_b}"

        guess = cm.block_overlap_length(block_b, block_a)
        assert (
            correct == guess
        ), f"Expected {correct} overlap but got {guess} for {block_b}, {block_a}"
