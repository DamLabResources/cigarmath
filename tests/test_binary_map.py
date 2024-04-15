import cigarmath as cm

def test_segments_to_binary():
    
    alns = [(30, cm.cigarstr2tup("30M"))]
    guess = cm.segments_to_binary(alns, max_genome_size = 100)
    assert sum(guess[:30]) == 0
    assert sum(guess[30:60]) == 30
    assert sum(guess[60:]) == 0
    
    alns = [(10, cm.cigarstr2tup("20M")),
            (50, cm.cigarstr2tup("20M"))
           ]
    
    guess = cm.segments_to_binary(alns, max_genome_size = 100)
    
    assert sum(guess[:10]) == 0
    assert sum(guess[10:30]) == 20
    
    assert sum(guess[30:50]) == 0
    assert sum(guess[50:70]) == 20
    
    assert sum(guess[70:]) == 0
    
    alns = [(10, cm.cigarstr2tup("20M10D5M")),
            (50, cm.cigarstr2tup("20M10D5M"))]
    guess = cm.segments_to_binary(alns, max_genome_size = 100, deletion_size=1)
    
    assert sum(guess[:10]) == 0
    assert sum(guess[10:30]) == 20
    assert sum(guess[30:40]) == 0
    assert sum(guess[40:45]) == 5
    
    assert sum(guess[45:50]) == 0
    assert sum(guess[50:70]) == 20
    assert sum(guess[70:80]) == 0
    assert sum(guess[80:85]) == 5
    
    assert sum(guess[85:]) == 0
    
    
def test_segments_to_binary_multi():
    
    alns1 = [(10, cm.cigarstr2tup("20M"))]
    alns2 = [(50, cm.cigarstr2tup("20M"))]
    
    # Process one alignment stream
    guess = cm.segments_to_binary(alns1, max_genome_size = 100)
    
    # Process a second alignment stream starting with the previous guess
    guess = cm.segments_to_binary(alns2, mapping = guess)
    
    assert sum(guess[:10]) == 0
    assert sum(guess[10:30]) == 20
    
    assert sum(guess[30:50]) == 0
    assert sum(guess[50:70]) == 20
    
    assert sum(guess[70:]) == 0