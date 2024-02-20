import cigarmath as cm
from cigarmath.defn import cigarstring_to_cigartuples as cigarstr2tup

def test_segments_to_binary():
    
    cigar, start = cigarstr2tup("30M"), 30
    guess = cm.plotting.segments_to_binary([cigar], [start], max_genome_size = 100)
    assert guess[:30].sum() == 0
    assert guess[30:60].sum() == 30
    assert guess[60:].sum() == 0
    
    
    cigars = [cigarstr2tup("20M")]*2
    starts = [10, 50]
    guess = cm.plotting.segments_to_binary(cigars, starts, max_genome_size = 100)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    
    assert guess[30:50].sum() == 0
    assert guess[50:70].sum() == 20
    
    assert guess[70:].sum() == 0
    
    
    cigars = [cigarstr2tup("20M10D5M")]*2
    starts = [10, 50]
    guess = cm.plotting.segments_to_binary(cigars, starts, max_genome_size = 100, deletion_size=1)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    assert guess[30:40].sum() == 0
    assert guess[40:45].sum() == 5
    
    assert guess[45:50].sum() == 0
    assert guess[50:70].sum() == 20
    assert guess[70:80].sum() == 0
    assert guess[80:85].sum() == 5
    
    assert guess[85:].sum() == 0