import numpy as np
import matplotlib.pyplot as plt

import cigarmath as cm

def test_segments_to_binary():
    
    alns = [(30, cm.cigarstr2tup("30M"))]
    guess = cm.segments_to_binary(alns, max_genome_size = 100)
    assert guess[:30].sum() == 0
    assert guess[30:60].sum() == 30
    assert guess[60:].sum() == 0
    
    alns = [(10, cm.cigarstr2tup("20M")),
            (50, cm.cigarstr2tup("20M"))
           ]
    
    guess = cm.segments_to_binary(alns, max_genome_size = 100)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    
    assert guess[30:50].sum() == 0
    assert guess[50:70].sum() == 20
    
    assert guess[70:].sum() == 0
    
    alns = [(10, cm.cigarstr2tup("20M10D5M")),
            (50, cm.cigarstr2tup("20M10D5M"))]
    guess = cm.segments_to_binary(alns, max_genome_size = 100, deletion_size=1)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    assert guess[30:40].sum() == 0
    assert guess[40:45].sum() == 5
    
    assert guess[45:50].sum() == 0
    assert guess[50:70].sum() == 20
    assert guess[70:80].sum() == 0
    assert guess[80:85].sum() == 5
    
    assert guess[85:].sum() == 0
    
    
def test_segments_to_binary_multi():
    
    alns1 = [(10, cm.cigarstr2tup("20M"))]
    alns2 = [(50, cm.cigarstr2tup("20M"))]
    
    # Process one alignment stream
    guess = cm.segments_to_binary(alns1, max_genome_size = 100)
    
    # Process a second alignment stream starting with the previous guess
    guess = cm.segments_to_binary(alns2, mapping = guess)
    
    assert guess[:10].sum() == 0
    assert guess[10:30].sum() == 20
    
    assert guess[30:50].sum() == 0
    assert guess[50:70].sum() == 20
    
    assert guess[70:].sum() == 0
    
    
def test_binary_map_plot():
    
    bmp = cm.plotting.BinaryMapPlot(np.random.rand(325, 1000)>0.2)
    
    map_ax = bmp.plot_map()
    hist_ax = bmp.plot_hist()
    
    assert isinstance(map_ax, plt.Axes)
    assert isinstance(hist_ax, plt.Axes)
    
def test_binary_map_plot_from_sam():
    
    stream = cm.io.segment_stream_pysam('tests/test_data/test.sam', 
                                    mode='r', progbar=True)
    
    bmp = cm.plotting.BinaryMapPlot.from_samstream(stream)
    
    map_ax = bmp.plot_map()
    hist_ax = bmp.plot_hist()
    
    assert isinstance(map_ax, plt.Axes)
    assert isinstance(hist_ax, plt.Axes)