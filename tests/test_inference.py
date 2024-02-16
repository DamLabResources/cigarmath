"""Test basic inference operations"""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from cigarmath.defn import cigarstring_to_cigartuples as cigarstr2tup

import cigarmath as cm

def check_cigartuples(guess, correct):
    "Compare two lists of cigartuples"
    assert len(guess) == len(correct)
    for (gop, gsz), (cop, csz) in zip(guess, correct):
        assert gop == cop
        assert gsz == csz