import pysam

import cigarmath as cm

def test_pysam_stream():
    
    stream = cm.io.segment_stream_pysam('tests/test_data/test.sam', 
                                        mode='r')
    
    for num, aln in enumerate(stream):
        assert isinstance(aln, pysam.AlignedSegment)
    
    assert num == 247 # Correct number of records
    
    
def test_pysam_stream_mapq():
    
    stream = cm.io.segment_stream_pysam('tests/test_data/test.sam', 
                                        mode='r',
                                        min_mapq=30)
    
    for num, aln in enumerate(stream):
        assert isinstance(aln, pysam.AlignedSegment)
    
    assert num == 198 # Correct number of records

def test_combined_segment_stream():

    stream = sorted(cm.io.segment_stream_pysam('tests/test_data/test.sam', mode='r'), key=lambda x: x.query_name)
    stream = cm.io.combined_segment_stream(stream)
    
    for num, (start, cigars, segments) in enumerate(stream):
        assert isinstance(segments[0], pysam.AlignedSegment)

    assert num == 192 # Correct number of records