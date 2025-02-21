# Cigar Math

A lightweight Python library for performing mathematical operations on CIGAR strings and tuples commonly used in genomic alignments.

## Overview

Cigar Math provides tools to handle and manipulate CIGAR strings/tuples - the compact representations of sequence alignments used in SAM/BAM files. It simplifies common operations like:

- Handling soft/hard clipping
- Finding overlapping alignment blocks
- Detecting deletion locations
- Converting between different alignment formats
- Mapping between reference and query coordinates

## Installation

```bash
pip install git+https://github.com/DamLabResources/cigarmath
```

## Basic Usage

```python
import cigarmath as cm

# Convert CIGAR string to tuples
cigartuples = cm.cigarstr2tup('3H4M1D3M2I3M4H')
reference_start = 3

# Get reference coordinates
ref_block = cm.reference_block(cigartuples, reference_start)
print(ref_block)  # (3, 14)

# Find mapping blocks (skipping deletions)
blocks = list(cm.reference_mapping_blocks(cigartuples, deletion_split=1))
print(blocks)  # [(0, 7), (10, 14), (20, 24)]

# Convert between coordinate spaces
query_pos = list(cm.reference2query(cigartuples, reference_start=2))
print(query_pos)  # [1, 2, 5, None, None, 6]
```

## Features

### Coordinate Mapping
- Map between reference and query coordinates
- Handle soft and hard clipping
- Find alignment blocks and gaps

### CIGAR Operations
- Parse and manipulate CIGAR strings/tuples
- Simplify extended CIGAR operations
- Collapse adjacent blocks

### Alignment Analysis
- Detect deletions and insertions
- Find overlapping regions
- Calculate alignment metrics

## Documentation

The `notebooks/` directory contains Jupyter notebooks with detailed examples:
- `block.ipynb`: Working with alignment blocks
- More examples coming soon...

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{cigar_math,
  author = {Dampier, Will and Klopfenstein, DV},
  title = {Cigar Math: A Python library for CIGAR string operations},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/DamLabResources/cigarmath}
}
```
