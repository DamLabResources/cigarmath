from cigarmath.defn import CONSUMES_REFERENCE
from cigarmath.defn import CONSUMES_QUERY

from cigarmath.clipping import left_clipping
from cigarmath.clipping import right_clipping
from cigarmath.clipping import declip



def cigartuples2pairs(cigartuples, reference_start = 0, verbose=False, clipping_fill = None):
    
    reference = []
    query = []
    
    left_clip = left_clipping(cigartuples)
    if left_clip:
        if verbose: print('Adding left clipping', left_clip)
        query += list(range(left_clip+1))
        reference = [clipping_fill]*left_clip
        if verbose: print('Currently', list(zip(query, reference)))
    
    for op, sz in declip(cigartuples):
        qstart = query[-1] if query else -1
        if verbose: print('Considering', op, sz)
        if op in CONSUMES_QUERY:
            if verbose: print('Consumes query, adding', list(range(qstart+1, qstart+sz+1)))
            query += list(range(qstart+1, qstart+sz+1))
        else:
            if verbose: print('Skips query, adding', [qstart]*sz)
            query += [query[-1]]*sz
        
        rstart = reference[-1] if (reference and reference[-1]) else reference_start-1
        if op in CONSUMES_REFERENCE:
            if verbose: print('Consumes Ref, adding', list(range(rstart+1, rstart+sz+1)))
            
            reference += list(range(rstart+1, rstart+sz+1))
        else:
            if verbose: print('Skips Ref, adding',  [rstart]*sz)
            reference += [rstart]*sz
        
        if verbose: print('Currently', list(zip(query, reference)))
            
    right_clip = right_clipping(cigartuples)
    if right_clip:
        if verbose: print('Adding right clipping:', right_clip)
        query += list(range(query[-1]+1, query[-1]+right_clip+1))
        reference += [clipping_fill]*right_clip
        if verbose: print('Currently', list(zip(query, reference)))
            
    return list(zip(query, reference))
    
    