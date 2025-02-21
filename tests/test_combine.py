"""Test combining CIGAR operations"""

from cigarmath.defn import cigarstr2tup, BAM_CDEL, BAM_CINS, BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CHARD_CLIP
from cigarmath.combine import (
    trim_alignment,
    combine_adjacent_alignments,
    combine_multiple_alignments
)

def test_trim_alignment():
    """Test trimming query bases from alignments"""
    
    # Test left trimming
    cigar = "5M1D3M"  # AAAAADCCC
    cigartups = cigarstr2tup(cigar)
    ref_start = 10
    
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=2)
    assert new_start == 12
    assert tuple(new_cigars) == ((0,3), (2,1), (0,3))
    
    # Test right trimming
    new_start, new_cigars = trim_alignment(ref_start, cigartups, right=2)
    assert new_start == 10
    assert tuple(new_cigars) == ((0,5), (2,1), (0,1))
    
    # Test both ends
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=2, right=2)
    assert new_start == 12
    assert tuple(new_cigars) == ((0,3), (2,1), (0,1))
    
    # Test with insertions
    cigar = "5M2I3M"  # AAAAAIICCC
    cigartups = cigarstr2tup(cigar)
    
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=6)
    assert new_start == 15
    assert tuple(new_cigars) == ((1,1),(0, 3))
    
    # Test no trimming
    new_start, new_cigars = trim_alignment(ref_start, cigartups)
    assert new_start == ref_start
    assert tuple(new_cigars) == tuple(cigartups)

    # Test spanning multiple blocks
    # RN : 0123456789---012345678901234567890
    # Q  : MMMMMMMMMMIIIMMM--MMMMMM
    # QN : 0123456789012345--678901

    cigar = '10M 3I 3M 2D 6M'
    ref_start = 0
    cigartups = cigarstr2tup(cigar)

    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=12)
    correct_cigars = [
            (BAM_CINS, 1),
            (BAM_CMATCH, 3),
            (BAM_CDEL, 2),
            (BAM_CMATCH, 6),
    ]
    assert new_start == 10
    assert tuple(new_cigars) == tuple(correct_cigars)
    
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=15)
    correct_cigars = [
            (BAM_CMATCH, 1),
            (BAM_CDEL, 2),
            (BAM_CMATCH, 6),
    ]
    assert new_start == 12
    assert tuple(new_cigars) == tuple(correct_cigars)

    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=18)
    correct_cigars = [
            (BAM_CMATCH, 4),
    ]
    print(new_cigars, correct_cigars)
    assert new_start == 17
    assert tuple(new_cigars) == tuple(correct_cigars)


def test_trim_alignment_with_clipping():
    """Test trimming query bases from alignments with clipping"""
    
    # Test adding soft clipping
    ref_start = 10
    cigartups = cigarstr2tup("5M1D3M")
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=2, right=2, add_clipping='soft')
    assert new_start == 12
    assert tuple(new_cigars) == ((BAM_CSOFT_CLIP, 2), (0,3), (2,1), (0,1), (BAM_CSOFT_CLIP, 2))
    
    # Test adding hard clipping
    new_start, new_cigars = trim_alignment(ref_start, cigartups, left=2, right=2, add_clipping='hard')
    assert new_start == 12
    assert tuple(new_cigars) == ((BAM_CHARD_CLIP, 2), (0,3), (2,1), (0,1), (BAM_CHARD_CLIP, 2))
    
    # Test invalid clipping type
    try:
        trim_alignment(ref_start, cigartups, left=2, add_clipping='invalid')
        assert False, "Should raise ValueError for invalid clipping type"
    except ValueError:
        pass


def test_combine_adjacent_alignments():
    """Test combining two adjacent alignments"""
    
    # REF: 0123456789012345678901234567890
    # Q1 :      AAAAASSS
    # Q2 :           SSSSSCCC
    # COR:      AAAAADDDDDCCC


    # Test with gap (deletion)
    first = (5, cigarstr2tup("5M3S"))
    second = (15, cigarstr2tup("5S3M"))

    correct_cigars = [
        (BAM_CMATCH, 5),
        (BAM_CDEL, 5),
        (BAM_CMATCH, 3),
    ]
    
    new_start, new_cigars = combine_adjacent_alignments(first, second)
    assert new_start == 5
    assert tuple(new_cigars) == tuple(correct_cigars)
    

    # Test with overlap
    # REF: 0123456789012345678|--|901234567890
    # Q1 :      AAAAAAAAAAAASSSSSS
    # Q2 :           CCC--CCCC|CC|CCCC
    # COR:      AAAAAAAAAAAACC|CC|CCCC


    # Test with gap (deletion)
    first = (5, cigarstr2tup("12M 17S"))
    second = (10, cigarstr2tup("12S 3M 2D 8M 2I 4M"))

    correct_cigars = [
        (BAM_CMATCH, 16),
        (BAM_CINS, 2),
        (BAM_CMATCH, 4),
    ]    
    
    new_start, new_cigars = combine_adjacent_alignments(first, second)
    assert new_start == 5
    assert tuple(new_cigars) == tuple(correct_cigars)

    # Test perfectly adjacent
    first = (10, cigarstr2tup("5M"))    # AAAAA
    second = (15, cigarstr2tup("3M"))   # CCC
    
    new_start, new_cigars = combine_adjacent_alignments(first, second)
    assert new_start == 10
    assert tuple(new_cigars) == ((0,8),)
    


def test_combine_multiple_alignments():
    """Test combining multiple alignments"""
    
    # Test basic combination

    alignments = [
        (5, cigarstr2tup("5M 2D 3M 18S")),  
        (20, cigarstr2tup("8S 4M 2I 2M 5S")),
        (29, cigarstr2tup("16S 5M 2I 3M"))
    ]

    # NUM: 012345678901234567890123--4567890123--4567890
    # Q1:       abcde--fgh
    # Q2:                      ijklMNop
    # Q3:                                 qrstUVwxyz
    # COR:      abcde--fgh-----ijklMNop---qrstUVwxyz

    new_start, new_cigars = combine_multiple_alignments(alignments)
    assert new_start == 5

    correct_cigars = [
        (BAM_CMATCH, 5),
        (BAM_CDEL, 2),
        (BAM_CMATCH, 3),

        (BAM_CDEL, 5),

        (BAM_CMATCH, 4),
        (BAM_CINS, 2),
        (BAM_CMATCH, 2),
        
        (BAM_CDEL, 3),

        (BAM_CMATCH, 5),
        (BAM_CINS, 2),
        (BAM_CMATCH, 3),
    ]

    assert tuple(new_cigars) == tuple(correct_cigars)
    
    overlapping = [
        (5, cigarstr2tup("5M 2D 3M 18S")),  
        (9, cigarstr2tup("8S 4M 2I 2M 10S")),
        (29, cigarstr2tup("16S 5M 2I 3M"))
    ]

    # NUM: 0123456789012345678901234567890123--4567890
    # Q1:       abcde--fgh
    # Q2:           ijklMNop
    # Q3:                               qrstuvwxyz
    # COR:      abcde--fghop------------qrstuvwxyz

    # Should fail with no overlap allowed
    try:
        combine_multiple_alignments(overlapping)
        assert False, "Should raise ValueError for overlap"
    except ValueError:
        pass
    

    correct_cigars = [
        (BAM_CMATCH, 5),
        (BAM_CDEL, 2),
        (BAM_CMATCH, 5),
        
        (BAM_CDEL, 12),

        (BAM_CMATCH, 5),
        (BAM_CINS, 2),
        (BAM_CMATCH, 3),
    ]

    # Should succeed with overlap allowed
    new_start, new_cigars = combine_multiple_alignments(overlapping, allowed_overlap=10)
    print(correct_cigars)
    print(new_cigars)
    
    assert new_start == 5
    assert tuple(new_cigars) == tuple(correct_cigars)
    
    # Test empty input
    try:
        combine_multiple_alignments([])
        assert False, "Should raise ValueError for empty input"
    except ValueError:
        pass
    
    # Test single alignment
    single = [(10, cigarstr2tup("5M"))]
    new_start, new_cigars = combine_multiple_alignments(single)
    assert new_start == 10
    assert tuple(new_cigars) == ((0,5),)
    
    # Test non-sequential reference positions
    non_sequential = [
        (20, cigarstr2tup("5M")),
        (10, cigarstr2tup("5M")),
    ]
    
    try:
        combine_multiple_alignments(non_sequential)
        assert False, "Should raise ValueError for non-sequential alignments"
    except ValueError:
        pass
