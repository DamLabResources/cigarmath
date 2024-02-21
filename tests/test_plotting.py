import cigarmath as cm

def test_segments_to_binary():
    
    alns = [(30, cm.cigarstr2tup("30M"))]
    guess = cm.plotting.segments_to_binary(alns, max_genome_size = 100)
    assert guess[:30].sum() == 0
    assert guess[30:60].sum() == 30
    assert guess[60:].sum() == 0
    
    alns = [(10, cm.cigarstr2tup("20M")),
            (50, cm.cigarstr2tup("20M"))
           ]
    
    guess = cm.plotting.segments_to_binary(alns, max_genome_size = 100)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    
    assert guess[30:50].sum() == 0
    assert guess[50:70].sum() == 20
    
    assert guess[70:].sum() == 0
    
    alns = [(10, cm.cigarstr2tup("20M10D5M")),
            (50, cm.cigarstr2tup("20M10D5M"))]
    guess = cm.plotting.segments_to_binary(alns, max_genome_size = 100, deletion_size=1)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    assert guess[30:40].sum() == 0
    assert guess[40:45].sum() == 5
    
    assert guess[45:50].sum() == 0
    assert guess[50:70].sum() == 20
    assert guess[70:80].sum() == 0
    assert guess[80:85].sum() == 5
    
    assert guess[85:].sum() == 0