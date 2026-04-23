CIGAR Format
============

CIGAR stands for Compact Idiosyncratic Gapped Alignment Report and is the
alignment edit string used in SAM/BAM files. A CIGAR encodes how a query/read
maps to a reference as a sequence of ``<length><operation>`` blocks.

Examples:

- ``10M``: 10 aligned positions (match/mismatch in SAM ``M`` semantics)
- ``4M1D3M2I3M``: aligned block, deletion, aligned block, insertion, aligned block
- ``3H4M1D3M2I3M4H``: hard clipping at both ends around the mapped region


Where CIGAR Is Used
-------------------

- SAM/BAM ``CIGAR`` field in read alignments.
- Coordinate conversion tools (query <-> reference position mapping).
- Coverage/pileup calculations and deletion/insertion-aware interval logic.
- Post-processing split alignments and trimming/merging alignments.

In this project, CIGAR operations are primarily represented as ``cigartuples``:

.. code-block:: python

    import cigarmath as cm

    cigartuples = cm.cigarstr2tup("3H4M1D3M2I3M4H")
    # [(5, 3), (0, 4), (2, 1), (0, 3), (1, 2), (0, 3), (5, 4)]

    cm.cigartup2str(cigartuples)
    # "3H4M1D3M2I3M4H"


Operation Semantics
-------------------

The core operations used by SAM/BAM and ``cigarmath``:

+-----------+------------------------------+----------------+---------------------+
| Op letter | Meaning                      | Consumes query | Consumes reference  |
+===========+==============================+================+=====================+
| ``M``     | Alignment match/mismatch     | Yes            | Yes                 |
+-----------+------------------------------+----------------+---------------------+
| ``I``     | Insertion to reference       | Yes            | No                  |
+-----------+------------------------------+----------------+---------------------+
| ``D``     | Deletion from reference      | No             | Yes                 |
+-----------+------------------------------+----------------+---------------------+
| ``N``     | Reference skip               | No             | Yes                 |
+-----------+------------------------------+----------------+---------------------+
| ``S``     | Soft clipping                | Yes            | No                  |
+-----------+------------------------------+----------------+---------------------+
| ``H``     | Hard clipping                | No             | No                  |
+-----------+------------------------------+----------------+---------------------+
| ``P``     | Padding                      | No             | No                  |
+-----------+------------------------------+----------------+---------------------+
| ``=``     | Sequence match               | Yes            | Yes                 |
+-----------+------------------------------+----------------+---------------------+
| ``X``     | Sequence mismatch            | Yes            | Yes                 |
+-----------+------------------------------+----------------+---------------------+

In ``cigarmath``, operation codes are integer indices into this operation list
(``M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8, B=9``), exposed in
``cigarmath.defn`` via constants like ``BAM_CMATCH``, ``BAM_CDEL``, and
``BAM_CSOFT_CLIP``.


Conventions in cigarmath
------------------------

- Coordinates are 0-based.
- Many coordinate mapping methods may return ``None`` at deleted or skipped positions.
- Clipping-aware operations are available (``left_clipping``, ``right_clipping``,
  ``declip``) to control whether clipped sequence contributes to downstream logic.
- ``M`` is treated as generic aligned positions (matches/mismatches), while
  extended semantics can use ``=`` and ``X``.


Canonical References
--------------------

- SAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- pysam ``cigartuples`` notes:
  https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

For practical examples using these semantics, continue to :doc:`usage`.
