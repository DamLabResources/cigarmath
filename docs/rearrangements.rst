Rearrangement inference
=======================

Split-read alignments from a single molecule can be classified by comparing
adjacent mapping segments along the original read. The decision tree below
summarizes how ``infer_rearrangements`` dispatches each segment pair.

.. image:: rearrangements.png
   :alt: Rearrangement inference decision tree
   :align: center

See :func:`cigarmath.rearrangement.infer_rearrangements` and
:func:`cigarmath.rearrangement.rearrangement_segment_stream`.

Diagram source
--------------

The canonical diagram source is ``rearrangements.dot`` (Graphviz). A Mermaid
version is kept in ``rearrangements.mmd`` for editors that prefer Mermaid syntax.

Regenerate the PNG after editing the decision tree::

    dot -Tpng docs/rearrangements.dot -o docs/rearrangements.png -Gdpi=150
