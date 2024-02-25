"""Various useful plotting operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

import numpy as np
import matplotlib.pyplot as plt

from cigarmath.block import reference_block
from cigarmath.block import reference_mapping_blocks


def segments_to_binary(alns, max_genome_size=10_000, deletion_size=50, mapping=None):
    """Given a list of segments create a binary array of all covered regions.
    
    POS0    000000000011111111112222222222333333333
    POS1    012345678901234567890123456789012345678
    REF     AAAAGACCCCCGACTAGCTAGCATGCTATCTAGCTAGCA
    QRY1    AAAAGA-----GAC      
    QRY2                           TGCTA---AGCTAG
    RES     111111000001110000000001111111111111100
    
    alns = [(0, '6M5D3M'),
            (23, '5M3D6M')]
    
    >> segments_to_binary(alns, max_genome_size=38, deletion_size=4)
    """
    
    if mapping is None:
        mapping = np.zeros(max_genome_size, dtype=bool)
    
    for start, cigar in alns:
        block_iter = reference_mapping_blocks(cigar,
                                              reference_start=start,
                                              deletion_split=deletion_size)
        for ref_start, ref_stop in block_iter:
            mapping[ref_start:ref_stop] = 1
        
    return mapping


class BinaryMapPlot:
    
    def __init__(self, binary_mat):
        self.binary_mat = binary_mat
    
    def __str__(self):
        return f'BinaryMapPlot(self.binary_mat = {self.binary_mat.shape})'
    
    def _repr__(self):
        return str(self)
    
    def plot_map(self, aspect = 'auto', cmap='binary', ax=None, **imshow_kws):
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
            
        
        ax.imshow(self.binary_mat, cmap=cmap, aspect=aspect, **imshow_kws)
        
        ax.set_ylabel('Reads')
        ax.set_xlabel('Genomic Position')
        
        return ax
    
    def plot_hist(self, bin_width=100, bins=None, ax=None, **hist_kws):
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
            
        if bins is None:
            bins = np.arange(0, self.binary_mat.shape[1], bin_width)
        
        ax.hist(self.binary_mat.sum(axis=1), bins = bins, **hist_kws)
        ax.set_ylabel('Reads')
        ax.set_xlabel('Mapped bases')
        
        return ax
    
    def dual_plot(self, bin_width=100, bins=None,  aspect = 'auto', cmap='binary', 
                  figsize = None,
                  hist_kws=None, imshow_kws=None):
        
        fig, (map_ax, hist_ax) = plt.subplots(1, 2, figsize=figsize,
                                              sharex=True, sharey=False)
        
        imshow_kws = imshow_kws if imshow_kws else {}
        self.plot_map(ax=map_ax, aspect = aspect, cmap=cmap, **imshow_kws)
        
        hist_kws = hist_kws if hist_kws else {}
        self.plot_hist(ax=hist_ax, bins=bins, bin_width=bin_width, **hist_kws)
        
        return fig, (map_ax, hist_ax)
    
    def save(self, handle):
        return np.save(handle, self.binary_mat)
    
    ## Constructors
    
    @staticmethod
    def load(handle):
        return BinaryMapPlot(np.load(handle))
    
    @staticmethod
    def from_samstream(stream, sort='mapped_bases',
                       max_genome_size=10_000,
                       deletion_size=50,
                       return_object=True):
        
        binary_matrix = BinaryMapPlot._stream2binary(stream, 
                                                     max_genome_size=max_genome_size,
                                                     deletion_size=deletion_size)
    
        if binary_matrix is None: return None
    
        if sort:
            binary_matrix = BinaryMapPlot._sort_matrix(binary_matrix, method = sort)
        
        return BinaryMapPlot(binary_matrix)
    
    
    ### HELPERS
    
    @staticmethod
    def _stream2binary(pysam_stream, max_genome_size=10_000, deletion_size=50):
        """Consumes a stream of pysam alingments and create a binary matrix."""
        
        binary_vectors = {}
        for segment in pysam_stream:
            if segment.query_name:
                qname = segment.query_name
                if segment.cigartuples:
                    binary_vectors[qname] = segments_to_binary([(segment.reference_start, segment.cigartuples)], 
                                                               max_genome_size=max_genome_size, 
                                                               deletion_size=deletion_size, 
                                                               mapping=binary_vectors.get(qname))
        
        if binary_vectors:
            return np.vstack(list(binary_vectors.values()))
    
    @staticmethod
    def _sort_matrix(matrix, method = 'mapped_bases'):
        
        if method == 'mapped_bases':
            return _sort_by_mapped_bases(matrix)
        raise ValueError(f'Did not understand sorting method: {method}')

        
        

def _sort_by_mapped_bases(matrix):
    "Sort the matrix by the number of aligned bases"
    
    mapped_bases = matrix.sum(axis=1)
    inds = np.argsort(mapped_bases)
    return matrix[inds, :]