=====
Usage
=====

This page focuses on practical workflows for working with CIGAR strings and
``cigartuples`` in Python.

A CIGAR string compactly describes how a query aligns to a reference (for
example ``3H4M1D3M2I3M4H``). In ``cigarmath``, most operations use
``cigartuples`` represented as ``(operation_code, length)`` pairs such as
``[(5, 3), (0, 4), (2, 1), ...]``. For a detailed explanation of the CIGAR
format and operation semantics, see :doc:`cigar`.

The examples below use top-level imports:

.. code-block:: python

    import cigarmath as cm


Quick Start
-----------

Parse a CIGAR string and get the mapped reference interval.

.. code-block:: python

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    reference_start = 3
    cm.reference_block(cigartuples, reference_start)

.. code-block:: text

    (3, 14)

Use this when you need a fast way to move from SAM-style CIGAR text to reference coordinates.


Coordinate Basics
-----------------

Get both reference and query intervals from the same alignment.

.. code-block:: python

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    ref_block = cm.reference_block(cigartuples, reference_start=3)
    qry_block = cm.query_block(cigartuples)
    qry_start = cm.query_start(cigartuples)
    ref_block, qry_block, qry_start

.. code-block:: text

    ((3, 14), (3, 15), 3)

Use this for coordinate bookkeeping before overlap, pileup, or interval conversion.


Clipping and Cleanup
--------------------

Inspect clipping and remove it before downstream operations.

.. code-block:: python

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    left = cm.left_clipping(cigartuples)
    right = cm.right_clipping(cigartuples)
    unclipped = cm.declip(cigartuples)
    hard = cm.is_hard_clipped(cigartuples)
    left, right, unclipped, hard

.. code-block:: text

    (3, 4, [(0, 4), (2, 1), (0, 3), (1, 2), (0, 3)], True)

Use this when clipping would otherwise shift query coordinates unexpectedly.


Mapping and Liftover
--------------------

Map positions between reference and query spaces.

.. code-block:: python

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    r2q = list(cm.reference2query(cigartuples, reference_start=3))
    q2r = list(cm.query2reference(cigartuples, reference_start=3))
    lift = list(cm.liftover(cigartuples, 3, 7, 8, 12, reference_start=3))
    r2q, q2r[:10], lift

.. code-block:: text

    ([3, 4, 5, 6, None, 7, 8, 9, 12, 13, 14],
     [None, None, None, 3, 4, 5, 6, 8, 9, 10],
     [3, 6, 7, 13])

`None` indicates unmapped positions (for example, deletions in one coordinate system).


Block Splitting and Overlap
---------------------------

Split alignments into mapped blocks around deletions and detect deletion intervals.

.. code-block:: python

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    mapped = list(cm.reference_mapping_blocks(cigartuples, reference_start=3, deletion_split=1))
    deletions = list(cm.reference_deletion_blocks(cigartuples, reference_start=3, min_size=1))
    overlap = cm.block_overlap_length((3, 14), (10, 20))
    mapped, deletions, overlap

.. code-block:: text

    ([(3, 7), (8, 14)], [(7, 8)], 4)

Use this to create interval-style representations from per-read alignments.


Combining and Trimming Alignments
---------------------------------

Trim query bases or merge neighboring alignments.

.. code-block:: python

    trimmed = cm.trim_alignment(
        10,
        ((0, 5), (2, 1), (0, 3)),
        left=2,
        right=1,
        add_clipping="soft",
    )

    combined = cm.combine_adjacent_alignments(
        (4, ((0, 4), (2, 1), (0, 3), (1, 2), (0, 3))),
        (12, ((0, 3), (0, 2))),
    )

    trimmed, combined

.. code-block:: text

    ((12, ((4, 2), (0, 3), (2, 1), (0, 2), (4, 1))),
     (4, ((0, 4), (2, 1), (0, 3), (1, 2), (0, 5))))

Use this when reconciling split or partially overlapping alignments.


Pileup and Base Counts
----------------------

Count per-position bases directly from CIGARs.

.. code-block:: python

    counts = cm.depth(
        [(0, 4), (2, 1), (0, 3), (1, 2), (0, 3)],
        query_sequence="AAAAACCGGCCC",
    )
    counts[0], counts[4], counts[9]

.. code-block:: text

    (Counter({'A': 1}), Counter({'-': 1}), Counter({'C': 1}))

Deletions are represented by `'-'` at reference positions.


Conversion Utilities
--------------------

Convert between CIGAR text, tuples, pair maps, and MSA-style strings.

.. code-block:: python

    cigar = cm.cigarstr2tup("1S2M2I1M2D1M2S")
    cigar_text = cm.cigartup2str(cigar)
    pairs = cm.cigartuples2pairs(cigar, reference_start=2)[:8]

    ref_msa = "AAAAGACCCCCGACTAGCTAGCATGCT----ATCTAGCTAGCA"
    qry_msa = "----AACCCCCGAC----TAGCATGCTTTTTATCTAGCT----"
    msa_start, msa_cigars = cm.msa2cigartuples(ref_msa, qry_msa)

    cigar_text, pairs, (msa_start, msa_cigars)

.. code-block:: text

    ('1S2M2I1M2D1M2S',
     [(0, None), (1, 2), (2, 3), (3, 3), (4, 3), (5, 4), (6, 5), (6, 6)],
     (4, [(0, 10), (2, 4), (0, 9), (1, 4), (0, 8)]))

Use these helpers when moving between aligner outputs and downstream analysis formats.


I/O Integration (Optional)
--------------------------

The `cigarmath.io` module can stream and combine alignments from SAM/BAM via `pysam`.
Install the optional dependency first; see :doc:`installation`.

Common entry points:

.. code-block:: python

    from cigarmath import io
    # io.segment_stream_pysam(...)
    # io.combined_segment_stream(...)


Common Pitfalls and Tips
------------------------

- Coordinates are 0-based throughout this library.
- `reference2query` and `query2reference` can return `None` at gaps/deletions.
- Run `declip(...)` when clipped ends are not part of your downstream coordinate model.
- For long deletions, tune `deletion_split` in `reference_mapping_blocks`.
- `combine_multiple_alignments` expects compatible query ordering and overlap behavior.
