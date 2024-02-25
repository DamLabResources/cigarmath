import cigarmath as cm

def test_aligned_pairs():
    
    cigar = '5M'
    pairs = cm.cigartuples2pairs(cm.cigarstr2tup(cigar))
    cor_pairs = [(0, 0),
                 (1, 1),
                 (2, 2),
                 (3, 3),
                 (4, 4),]
    assert list(pairs) == cor_pairs
    
    cigar = '5M'
    pairs = cm.cigartuples2pairs(cm.cigarstr2tup(cigar), reference_start=5)
    cor_pairs = [(0, 5),
                 (1, 6),
                 (2, 7),
                 (3, 8),
                 (4, 9),]
    assert list(pairs) == cor_pairs
    
    # REF 0123456--789
    # QRY      0123456
    
    cigar = '2M2I3M'
    pairs = cm.cigartuples2pairs(cm.cigarstr2tup(cigar), reference_start=5)
    cor_pairs = [(0, 5),
                 (1, 6),
                 (2, 6),
                 (3, 6),
                 (4, 7),
                 (5, 8),
                 (6, 9),
                ]
    assert list(pairs) == cor_pairs
    
    # REF 0123456789
    # QRY  01--234
    
    cigar = '2M2D3M'
    pairs = cm.cigartuples2pairs(cm.cigarstr2tup(cigar), reference_start=1)
    cor_pairs = [(0, 1),
                 (1, 2),
                 (1, 3),
                 (1, 4),
                 (2, 5),
                 (3, 6),
                 (4, 7),
                ]
    assert list(pairs) == cor_pairs

    cigar = '2S5M3H'
    pairs = cm.cigartuples2pairs(cm.cigarstr2tup(cigar), reference_start=5, verbose=True)
    cor_pairs = [(0, None),
                 (1, None),
                 (2, 5),
                 (3, 6),
                 (4, 7),
                 (5, 8),
                 (6, 9),
                 (7, None),
                 (8, None),
                 (9, None),
                ]
    assert list(pairs) == cor_pairs