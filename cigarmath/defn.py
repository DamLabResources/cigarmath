"""CIGAR Strings (Compact Idiosyncratic Gapped Alignment Report"""
# https://samtools.github.io/hts-specs/SAMv1.pdf

from collections import namedtuple

__copyright__ = (
    "Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved"
)
__author__ = "DV Klopfenstein, PhD"

NTO = namedtuple("CigarLetter", "Op desc consumes_query consumes_ref")

#                                                                            consumes
#              Op    desc                                                   query  ref
#              ---   -----------------------------------------------------  ------ -----
NTS = [
    NTO._make(
        ["M", "alignment match (can be a sequence match or mismatch)", True, True]
    ),
    NTO._make(["I", "insertion to the reference", True, False]),
    NTO._make(["D", "deletion form the reference", False, True]),
    NTO._make(["N", "skipped region from the reference", False, True]),
    NTO._make(["S", "soft clipping (clipped sequences present in SEQ)", True, False]),
    NTO._make(
        ["H", "hard clipping (clipped sequences NOT present in SEQ)", False, False]
    ),
    NTO._make(
        ["P", "padding (silent deletion from the padded reference)", False, False]
    ),
    NTO._make(["=", "sequnce match", True, True]),
    NTO._make(["X", "sequence mismatch", True, True]),
    # https://sourceforge.net/p/samtools/mailman/message/29373646/
    NTO._make(["B", 'proposed "backwards" skip operator, like N', False, False]),
    # https://samtools.github.io/hts-specs/SAMtags.pdf
    # Number of differences ([mismatches plus inserted and deleted bases]) betweet the
    # sequence and the reference, counting only A, C, G, and T bases in seq/ref
    # as potential matches, with everything else being a mismatch.
    # NOTE: Historically ill-defined w/both data and tools disagreeing with defn.
    NTO._make(["T", "NM tag: Edit distance to the reference", False, False]),
]

NtoRQC = namedtuple("CigRQC", "Op consumes_ref consumes_query")
NTS_RQC = [NtoRQC._make([nt.Op, nt.consumes_ref, nt.consumes_query]) for nt in NTS]


def get_ntrqc_qty(cigartuples):
    """Convert BAM AlignedSequence cigartuples to RefQryCIGAR nt2cnt"""
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    #     "alignment is returned as a list of tuples of (operation, length)"
    return [(NTS_RQC[i], n) for i, n in cigartuples]


CIGAR_HDRS = [nt.Op for nt in NTS]
CIGAR_SET = set(CIGAR_HDRS)
CIGAR_VALID_CHRS = set(str(n) for n in range(10)).union(CIGAR_SET)

CIGAR2BAM = {nt.Op: i for i, nt in enumerate(NTS)}

CONSUMES_REFERENCE = {
    bam_num for bam_num, cigar_letter in enumerate(NTS) if cigar_letter.consumes_ref
}
CONSUMES_QUERY = {
    bam_num for bam_num, cigar_letter in enumerate(NTS) if cigar_letter.consumes_query
}


def cigarstr2tup(cigarstring):
    """Create cigartuples from cigarstring"""
    if cigarstring is None:
        return None
    cigartuples = []
    pta = 0
    ##print(f'CIGAR({cigarstring})')
    for idx, letter in enumerate(cigarstring):
        if letter in CIGAR_SET:
            ##print(f'LETTER({letter}) VAL({cigarstring[pta:idx]})')
            cigartuples.append((CIGAR2BAM[letter], int(cigarstring[pta:idx])))
            pta = idx + 1
    return cigartuples


def cigartup2str(cigartuples, sep=""):
    """Get cigarstring from cigartuples"""
    return sep.join(f"{cnt}{CIGAR_HDRS[op]}" for op, cnt in cigartuples)


BAM_CMATCH = 0  # M
BAM_CINS = 1  # I
BAM_CDEL = 2  # D
BAM_CREF_SKIP = 3  # N
BAM_CSOFT_CLIP = 4  # S
BAM_CHARD_CLIP = 5  # H
BAM_CPAD = 6  # P
BAM_CEQUAL = 7  # =
BAM_CDIFF = 8  # X
BAM_CBACK = 9  # B


# Copyright (C) 2022-present, Dampier & DV Klopfenstein, PhD. All rights reserved
