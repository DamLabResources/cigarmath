"""Infer chromosomal rearrangements from split-read alignments."""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from collections import namedtuple
from dataclasses import dataclass
from itertools import groupby
from typing import Dict, Iterator, List, Optional, Tuple, TYPE_CHECKING

from cigarmath.block import query_block, reference_block
from cigarmath.defn import (
    BAM_CDEL,
    BAM_CINS,
    BAM_CREF_SKIP,
    CigarTuples,
    CONSUMES_QUERY,
    CONSUMES_REFERENCE,
)
from cigarmath.inference import inferred_query_sequence_length

if TYPE_CHECKING:
    import pysam

Mapping = Tuple[str, int, bool, CigarTuples]

RearrangementEvent = namedtuple(
    "RearrangementEvent",
    [
        "event_type",
        "segment_indices",
        "left_reference",
        "right_reference",
        "reference_size",
        "query_size",
        "strands",
    ],
)


def _rearrangement_event_str(event: RearrangementEvent) -> str:
    left_name, left_pos = event.left_reference
    right_name, right_pos = event.right_reference
    if event.event_type == "TRA":
        ref_part = "ref=N/A"
    elif event.event_type == "INTRA_INS":
        ref_part = f"ref={left_pos}"
        return (
            f"{event.event_type} {left_name}:{left_pos} @ q[{right_pos}) "
            f"(query=+{event.query_size})"
        )
    elif event.event_type == "INTRA_DEL":
        ref_part = f"ref +{event.reference_size}"
        return (
            f"{event.event_type} {left_name}:{left_pos}-{right_pos} @ q[{event.query_size}) "
            f"({ref_part})"
        )
    else:
        ref_part = f"ref={event.reference_size:+d}"
    return (
        f"{event.event_type} {left_name}:{left_pos}->{right_name}:{right_pos} "
        f"({ref_part}, query=+{event.query_size})"
    )


RearrangementEvent.__str__ = _rearrangement_event_str  # type: ignore[method-assign]


@dataclass
class _Segment:
    index: int
    ref_name: str
    ref_start: int
    ref_end: int
    is_reverse: bool
    q_start: int
    q_end: int


def _read_length_from_mappings(mappings: List[Mapping]) -> int:
    """Infer full read length as the maximum IQL across split mappings.

    One supplementary alignment may cover only part of the read; the longest
    implied query length is treated as the molecule length for coordinate projection.

    S0  IQL=1000   S1  IQL=1000   -> read_len=1000
    S0  IQL=800    S1  IQL=1000   -> read_len=1000
    """
    if not mappings:
        return 0
    return max(
        inferred_query_sequence_length(cigartuples)
        for _, _, _, cigartuples in mappings
    )


def _normalize_segments(
    mappings: List[Mapping],
    min_segment_length: int,
) -> Tuple[int, List[_Segment]]:
    """Project each mapping to read-ordered segments with reference and query spans.

    Mappings are sorted by projected query start. Reverse-strand records are flipped
    into read order (5'->3' along the sequenced molecule). Short alignments below
    ``min_segment_length`` on the reference are dropped.

    READ  ----[====S0====]----[==S1==]----
    REF       chr1:1000-1400    chr1:5000-5200

    Returns (read_len, segments).
    """
    if not mappings:
        return 0, []

    read_len = _read_length_from_mappings(mappings)

    segments: List[_Segment] = []
    for index, (ref_name, ref_start, is_reverse, cigartuples) in enumerate(mappings):
        ref_start_pos, ref_end_pos = reference_block(cigartuples, ref_start)
        if ref_end_pos - ref_start_pos < min_segment_length:
            continue

        q_start_seg, q_end_seg = query_block(cigartuples)
        if is_reverse:
            q_start = read_len - q_end_seg
            q_end = read_len - q_start_seg
        else:
            q_start, q_end = q_start_seg, q_end_seg

        segments.append(
            _Segment(
                index=index,
                ref_name=ref_name,
                ref_start=ref_start_pos,
                ref_end=ref_end_pos,
                is_reverse=is_reverse,
                q_start=q_start,
                q_end=q_end,
            )
        )

    segments.sort(key=lambda s: s.q_start)
    return read_len, segments


def _strand_symbol(is_reverse: bool) -> str:
    return "-" if is_reverse else "+"


def _segment_line(
    segment: _Segment,
    event_suffix: str = "",
    q_start: Optional[int] = None,
    q_end: Optional[int] = None,
    ref_start: Optional[int] = None,
    ref_end: Optional[int] = None,
) -> str:
    q_lo = segment.q_start if q_start is None else q_start
    q_hi = segment.q_end if q_end is None else q_end
    ref_lo = segment.ref_start if ref_start is None else ref_start
    ref_hi = segment.ref_end if ref_end is None else ref_end
    line = (
        f"[S{segment.index} {_strand_symbol(segment.is_reverse)}] "
        f"q[{q_lo},{q_hi}) "
        f"ref:{segment.ref_name}:{ref_lo}-{ref_hi}"
    )
    if event_suffix:
        line = f"{line} {event_suffix}"
    return line


def _intra_event_read_anchor(event: RearrangementEvent) -> int:
    if event.event_type == "INTRA_INS":
        return event.right_reference[1]
    return event.query_size


def _append_segment_with_intra_events(
    lines: List[str],
    segment: _Segment,
    intra_events: List[RearrangementEvent],
    event_suffix: str = "",
) -> None:
    if not intra_events:
        lines.append(_segment_line(segment, event_suffix))
        return

    sorted_events = sorted(intra_events, key=_intra_event_read_anchor)
    remaining_q = segment.q_start
    remaining_ref = segment.ref_start
    suffix_on_last = event_suffix

    for event in sorted_events:
        if event.event_type == "INTRA_DEL":
            q_at = event.query_size
            ref_del_start = event.left_reference[1]
            ref_del_end = event.right_reference[1]
            if q_at > remaining_q:
                lines.append(
                    _segment_line(
                        segment,
                        q_start=remaining_q,
                        q_end=q_at,
                        ref_start=remaining_ref,
                        ref_end=ref_del_start,
                    )
                )
            lines.append(_intra_segment_line(event))
            remaining_q = q_at
            remaining_ref = ref_del_end
        elif event.event_type == "INTRA_INS":
            q_at = event.right_reference[1]
            ref_at = event.left_reference[1]
            if q_at > remaining_q:
                lines.append(
                    _segment_line(
                        segment,
                        q_start=remaining_q,
                        q_end=q_at,
                        ref_start=remaining_ref,
                        ref_end=ref_at,
                    )
                )
            lines.append(_intra_segment_line(event))
            remaining_q = q_at + event.query_size
            remaining_ref = ref_at

    if remaining_q < segment.q_end or remaining_ref < segment.ref_end:
        lines.append(
            _segment_line(
                segment,
                event_suffix=suffix_on_last,
                q_start=remaining_q,
                q_end=segment.q_end,
                ref_start=remaining_ref,
                ref_end=segment.ref_end,
            )
        )
    elif suffix_on_last:
        lines.append(_segment_line(segment, event_suffix=suffix_on_last))


def _event_suffix(event: RearrangementEvent) -> str:
    if event.event_type == "TRA":
        if event.query_size:
            return f"| TRA query +{event.query_size}"
        return "| TRA"
    return f"| {event.event_type} ref {event.reference_size:+d}"


def _between_segment_line(event: RearrangementEvent) -> str:
    if event.event_type == "INS":
        return f"[INS] query +{event.query_size} bp"
    if event.event_type == "DEL":
        return f"[DEL] ref +{event.reference_size} bp"
    raise ValueError(f"unexpected between-segment event type: {event.event_type}")


def _intra_segment_line(event: RearrangementEvent) -> str:
    if event.event_type == "INTRA_INS":
        return f"[INTRA_INS] query +{event.query_size} bp"
    if event.event_type == "INTRA_DEL":
        return f"[INTRA_DEL] ref +{event.reference_size} bp"
    raise ValueError(f"unexpected intra-segment event type: {event.event_type}")


def _read_query_coord(raw_query_pos: int, is_reverse: bool, read_len: int) -> int:
    """Map a CIGAR-walk query coordinate into read order (5'->3' along the read).

    Forward alignments use raw positions unchanged. Reverse alignments mirror
    against ``read_len`` so intra-segment anchors match ``_normalize_segments``.

    read_len=1000, reverse, raw_q=700  ->  read_q=300
    """
    if is_reverse:
        return read_len - raw_query_pos
    return raw_query_pos


def _make_intra_event(
    event_type: str,
    segment_index: int,
    ref_name: str,
    ref_start: int,
    ref_end: int,
    read_query_anchor: int,
    reference_size: int,
    query_size: int,
    is_reverse: bool,
) -> RearrangementEvent:
    """Build an intra-segment event tied to a single mapping index.

    ``INTRA_DEL``: ``reference_size`` is deletion length; ``query_size`` is the
    read-coordinate anchor before the deleted reference bases.

    ``INTRA_INS``: ``query_size`` is insertion length; ``right_reference[1]`` stores
    the read-coordinate anchor; ``left_reference[1]`` is the reference position.
    """
    if event_type == "INTRA_INS":
        return RearrangementEvent(
            event_type=event_type,
            segment_indices=(segment_index, segment_index),
            left_reference=(ref_name, ref_start),
            right_reference=(ref_name, read_query_anchor),
            reference_size=reference_size,
            query_size=query_size,
            strands=(is_reverse, is_reverse),
        )
    return RearrangementEvent(
        event_type=event_type,
        segment_indices=(segment_index, segment_index),
        left_reference=(ref_name, ref_start),
        right_reference=(ref_name, ref_end),
        reference_size=reference_size,
        query_size=read_query_anchor,
        strands=(is_reverse, is_reverse),
    )


def _detect_intra_segment_indels(
    mappings: List[Mapping],
    min_event_size: int,
) -> List[RearrangementEvent]:
    """Detect large ``I``/``D`` (and ``N``) operations inside a single alignment CIGAR.

    These are indels within one SAM record, not gaps between supplementary segments.
    Emits ``INTRA_INS`` / ``INTRA_DEL`` with ``segment_indices=(i, i)``.

    READ  ----AAAA----IIII----CCCC----
    REF       AAAAGGGG----CCCC
    CGS       4M   200I   4M        -> INTRA_INS +200 at query anchor after first 4M
    CGS       4M   500D   4M        -> INTRA_DEL +500 at matching ref coordinates
    """
    if not mappings:
        return []

    read_len = _read_length_from_mappings(mappings)
    del_ops = {BAM_CDEL, BAM_CREF_SKIP}
    events: List[RearrangementEvent] = []

    for index, (ref_name, ref_start, is_reverse, cigartuples) in enumerate(mappings):
        ref_pos = ref_start
        query_pos = 0
        for op, size in cigartuples:
            if op == BAM_CINS and size >= min_event_size:
                read_q = _read_query_coord(query_pos, is_reverse, read_len)
                events.append(
                    _make_intra_event(
                        "INTRA_INS",
                        index,
                        ref_name,
                        ref_pos,
                        ref_pos,
                        read_q,
                        0,
                        size,
                        is_reverse,
                    )
                )
                query_pos += size
            elif op in del_ops and size >= min_event_size:
                read_q = _read_query_coord(query_pos, is_reverse, read_len)
                events.append(
                    _make_intra_event(
                        "INTRA_DEL",
                        index,
                        ref_name,
                        ref_pos,
                        ref_pos + size,
                        read_q,
                        size,
                        0,
                        is_reverse,
                    )
                )
                ref_pos += size
            else:
                if op in CONSUMES_REFERENCE:
                    ref_pos += size
                if op in CONSUMES_QUERY:
                    query_pos += size

    events.sort(
        key=lambda event: (
            event.segment_indices[0],
            event.right_reference[1]
            if event.event_type == "INTRA_INS"
            else event.query_size,
        )
    )
    return events


def format_read_rearrangement_summary(
    query_name: str,
    mappings: List[Mapping],
    events: List[RearrangementEvent],
    show_clip_threshold: int = 30,
    min_segment_length: int = 1,
) -> str:
    """ASCII bracket summary of split-read segments and inferred rearrangements."""
    read_len, segments = _normalize_segments(mappings, min_segment_length)
    lines = [
        f"read {query_name} len={read_len} segments={len(segments)} events={len(events)}"
    ]

    between_by_pair = {
        event.segment_indices: event
        for event in events
        if event.event_type in ("INS", "DEL")
    }
    intra_by_segment: dict[int, List[RearrangementEvent]] = {}
    for event in events:
        if event.event_type not in ("INTRA_INS", "INTRA_DEL"):
            continue
        intra_by_segment.setdefault(event.segment_indices[0], []).append(event)

    downstream_labels: dict[int, List[str]] = {}
    for event in events:
        if event.event_type in ("INS", "DEL", "INTRA_INS", "INTRA_DEL"):
            continue
        _, right_index = event.segment_indices
        downstream_labels.setdefault(right_index, []).append(_event_suffix(event))

    if segments and segments[0].q_start > show_clip_threshold:
        lines.append(f"[Clip] q[0, {segments[0].q_start})")

    for idx, segment in enumerate(segments):
        suffix = ""
        labels = downstream_labels.get(segment.index)
        if labels:
            suffix = " ".join(labels)
        _append_segment_with_intra_events(
            lines,
            segment,
            intra_by_segment.get(segment.index, []),
            event_suffix=suffix,
        )

        if idx + 1 < len(segments):
            next_segment = segments[idx + 1]
            pair = (segment.index, next_segment.index)
            between_event = between_by_pair.get(pair)
            if between_event is not None:
                lines.append(_between_segment_line(between_event))

    if segments and read_len - segments[-1].q_end > show_clip_threshold:
        lines.append(f"[Clip] q[{segments[-1].q_end}, {read_len})")

    return "\n".join(lines)


def _reference_proximity(left: _Segment, right: _Segment) -> int:
    """Minimum distance between reference endpoints of two adjacent segments.

    Used for inversion breakpoint clustering: small values mean the splits likely
    bracket one rearrangement junction.

    REF  ..........[====left====]....[====right====].......
                      ^left_end   right_start^
    proximity = min(|right_start - left_end|, |left_start - right_end|)

    left ends 1400, right starts 1380  -> proximity=20
    """
    return min(
        abs(right.ref_start - left.ref_end),
        abs(left.ref_start - right.ref_end),
    )


def _reference_overlap(left: _Segment, right: _Segment) -> int:
    """Length of reference interval shared by two segments (0 if none).

    REF  ....[====left====]
    REF  ........[====right====]...
    overlap = min(left_end, right_end) - max(left_start, right_start)

    left 1000-1400, right 1280-1580  -> overlap=120
    """
    return min(left.ref_end, right.ref_end) - max(left.ref_start, right.ref_start)


def _make_event(
    event_type: str,
    left: _Segment,
    right: _Segment,
    reference_size: Optional[int],
    query_size: int,
) -> RearrangementEvent:
    """Build a between-segment ``RearrangementEvent`` from an adjacent segment pair.

    ``left`` is upstream on the read; ``right`` is downstream. Breakpoints are stored
    as ``left_reference=(ref, left.ref_end)`` and ``right_reference=(ref, right.ref_start)``.
    """
    return RearrangementEvent(
        event_type=event_type,
        segment_indices=(left.index, right.index),
        left_reference=(left.ref_name, left.ref_end),
        right_reference=(right.ref_name, right.ref_start),
        reference_size=reference_size,
        query_size=query_size,
        strands=(left.is_reverse, right.is_reverse),
    )


def _detect_translocation(left: _Segment, right: _Segment) -> RearrangementEvent:
    """Classify adjacent segments on different contigs as translocation (TRA).

    READ  ----[====S0====]----[====S1====]----
    REF   chr1:1000-1400       chr7:80000-80400

    Always returns TRA; ``reference_size`` is None and ``query_size`` is the
    intervening query gap (if any).
    """
    query_gap = max(0, right.q_start - left.q_end)
    return _make_event("TRA", left, right, None, query_gap)


def _detect_inversion(
    left: _Segment,
    right: _Segment,
    max_breakpoint_distance: int,
    min_embedded_inversion_size: int = 150,
    min_embedded_overlap_fraction: float = 0.5,
    max_junction_slop: int = 20,
) -> Optional[RearrangementEvent]:
    """Detect inversion (INV) when adjacent segments map to opposite strands.

    Path 1 — close breakpoints (classic split-read inversion):

    READ  ----[====S0 +====]----[====S1 -====]----
    REF        chr1:1000-1400     chr1:1380-1780
    proximity = |1380 - 1400| = 20  <= max_breakpoint_distance  -> INV

    Path 2 — embedded inversion (overlapping +/- blocks, distant endpoints):

    READ  ----[========S0 +========]----
    READ              [====S1 -====]
    REF   chr1:2100-5500 (forward span contains reverse block)
    large ref overlap, query junction slop <= max_junction_slop  -> INV

    Returns None if strands match or neither path qualifies.
    """
    query_junction = right.q_start - left.q_end
    query_gap = max(0, query_junction)

    if _reference_proximity(left, right) <= max_breakpoint_distance:
        return _make_event(
            "INV", left, right, right.ref_start - left.ref_end, query_gap
        )

    ref_overlap = _reference_overlap(left, right)
    if ref_overlap < min_embedded_inversion_size:
        return None
    left_span = left.ref_end - left.ref_start
    right_span = right.ref_end - right.ref_start
    overlap_fraction = ref_overlap / max(1, min(left_span, right_span))
    if overlap_fraction < min_embedded_overlap_fraction:
        return None
    if abs(query_junction) > max_junction_slop:
        return None
    return _make_event("INV", left, right, -ref_overlap, query_gap)


def reference_lengths_from_pysam_header(header: dict) -> Dict[str, int]:
    """Build contig name -> length from a pysam header dict."""
    return {entry["SN"]: entry["LN"] for entry in header.get("SQ", [])}


def _is_genomic_reference_wrap(
    left: _Segment,
    right: _Segment,
    reference_size: int,
    segments: List[_Segment],
    reference_lengths: Optional[Dict[str, int]],
    end_frac_threshold: float = 0.75,
    start_frac_threshold: float = 0.25,
    min_overlap_genome_fraction: float = 0.5,
    min_read_span_genome_fraction: float = 0.7,
) -> bool:
    """Return True when a negative reference jump is linear-genome wrap, not tandem DUP.

    Same-strand split where the downstream segment maps far upstream on the reference
    (HIV 3' body then 5' terminus on HXB2F), not a local duplicate.

    READ  ----[======== late genome S0 ========][S1 proximal]----
    REF   HXB2F:2100-8582  then  HXB2F:54-180   (reference_size ~ -8500)

    Signals (any can qualify):
    - |reference_size| larger than both segment spans
    - left ends in distal fraction of chromosome, right starts in proximal fraction
    - read spans most of the chromosome with both termini represented

    >>> _is_genomic_reference_wrap(...)  # True for HXB2F LTR-wrap-like pairs
    """
    left_span = left.ref_end - left.ref_start
    right_span = right.ref_end - right.ref_start
    if abs(reference_size) > max(left_span, right_span):
        return True

    ref_len = None if reference_lengths is None else reference_lengths.get(left.ref_name)
    if ref_len is None or ref_len <= 0:
        return False

    left_frac = left.ref_end / ref_len
    right_frac = right.ref_start / ref_len
    overlap_frac = abs(reference_size) / ref_len
    if (
        left_frac >= end_frac_threshold
        and right_frac <= start_frac_threshold
        and overlap_frac >= min_overlap_genome_fraction
    ):
        return True

    same_chr = [segment for segment in segments if segment.ref_name == left.ref_name]
    if len(same_chr) < 2:
        return False
    span_start = min(segment.ref_start for segment in same_chr)
    span_end = max(segment.ref_end for segment in same_chr)
    if (span_end - span_start) / ref_len <= min_read_span_genome_fraction:
        return False
    has_proximal = any(
        segment.ref_start < ref_len * start_frac_threshold for segment in same_chr
    )
    has_distal = any(
        segment.ref_end > ref_len * end_frac_threshold for segment in same_chr
    )
    return has_proximal and has_distal


def _detect_duplication(
    left: _Segment,
    right: _Segment,
    min_event_size: int,
    segments: List[_Segment],
    reference_lengths: Optional[Dict[str, int]],
) -> Optional[RearrangementEvent]:
    """Detect tandem duplication (DUP) or reference wrap (REF_WRAP) on the same strand.

    Requires the downstream segment to start well before the upstream segment ends
    on the reference (negative ``reference_size``).

    READ  ----[====S0====]----[====S1====]----
    REF   chr1:1000-2000        chr1:1800-2800
    reference_size = 1800 - 2000 = -200  -> DUP (local overlap)

    READ  ----[======== S0 ========][S1 5' terminus]----
    REF   HXB2F:2100-8582  ->  HXB2F:54-180
    reference_size ~ -8500, genome-wrap heuristics  -> REF_WRAP

    Returns None if ``reference_size >= 0`` or event is smaller than ``min_event_size``.
    """
    reference_size = right.ref_start - left.ref_end
    if reference_size >= 0:
        return None
    if abs(reference_size) < min_event_size:
        return None
    query_gap = max(0, right.q_start - left.q_end)
    if _is_genomic_reference_wrap(
        left, right, reference_size, segments, reference_lengths
    ):
        return _make_event("REF_WRAP", left, right, reference_size, query_gap)
    return _make_event("DUP", left, right, reference_size, query_gap)


def _detect_deletion(
    left: _Segment,
    right: _Segment,
    min_event_size: int,
) -> Optional[RearrangementEvent]:
    """Detect deletion (DEL) between colinear same-strand segments.

    Reference coordinates jump forward farther than the read advances: sequence
    missing on the reference between split alignments, not within either CIGAR.

    READ  ----[====S0====]----[====S1====]----
    REF   chr1:1000-1500  ....gap....  chr1:11500-12000
    ref_gap=10000, query_gap=0  -> DEL +10000 bp

    Returns None if the reference gap is smaller than ``min_event_size`` or if the
    query gap is as large as the reference gap (gaps dominated by inserted sequence).
    """
    ref_gap = right.ref_start - left.ref_end
    query_gap = max(0, right.q_start - left.q_end)
    if ref_gap < min_event_size or query_gap >= ref_gap:
        return None
    return _make_event("DEL", left, right, ref_gap, query_gap)


def _detect_insertion(
    left: _Segment,
    right: _Segment,
    min_event_size: int,
) -> Optional[RearrangementEvent]:
    """Detect insertion (INS) between colinear same-strand segments.

    The read advances farther than reference coordinates: sequence present in the
    molecule between split alignments but not accounted for by reference span.

    READ  ----[====S0====]--extra--[====S1====]----
    REF   chr1:1000-1500              chr1:1500-1900
    query_gap=100, ref_gap=0  -> INS +100 bp

    Returns None if the query gap is smaller than ``min_event_size`` or if the
    reference gap is at least as large as the query gap (deletion-like spacing).
    """
    ref_gap = right.ref_start - left.ref_end
    query_gap = max(0, right.q_start - left.q_end)
    if query_gap < min_event_size or ref_gap >= query_gap:
        return None
    return _make_event("INS", left, right, ref_gap, query_gap)


def infer_rearrangements(
    mappings: List[Mapping],
    max_breakpoint_distance: int = 1000,
    min_segment_length: int = 50,
    min_event_size: int = 30,
    reference_lengths: Optional[Dict[str, int]] = None,
    min_embedded_inversion_size: int = 150,
    min_embedded_overlap_fraction: float = 0.5,
    max_junction_slop: int = 20,
) -> List[RearrangementEvent]:
    """Classify CIGAR-internal and between-segment rearrangements for one read.

    First pass: ``INTRA_INS`` / ``INTRA_DEL`` from large ``I`` and ``D`` ops in each
    mapping CIGAR. Second pass: compare adjacent query-ordered segments for
    TRA, INV, DUP/REF_WRAP, DEL, and INS.

    READ  ----[S0]----[S1]----   (two supplementary alignments, one molecule)
    """
    _, segments = _normalize_segments(mappings, min_segment_length)
    events: List[RearrangementEvent] = _detect_intra_segment_indels(
        mappings, min_event_size
    )

    if len(segments) < 2:
        return events

    for left, right in zip(segments, segments[1:]):
        if left.ref_name != right.ref_name:
            events.append(_detect_translocation(left, right))
        elif left.is_reverse != right.is_reverse:
            event = _detect_inversion(
                left,
                right,
                max_breakpoint_distance,
                min_embedded_inversion_size=min_embedded_inversion_size,
                min_embedded_overlap_fraction=min_embedded_overlap_fraction,
                max_junction_slop=max_junction_slop,
            )
            if event is not None:
                events.append(event)
        elif right.ref_start < left.ref_end - max_breakpoint_distance:
            event = _detect_duplication(
                left,
                right,
                min_event_size,
                segments,
                reference_lengths,
            )
            if event is not None:
                events.append(event)
        else:
            ref_gap = right.ref_start - left.ref_end
            query_gap = max(0, right.q_start - left.q_end)
            if ref_gap >= min_event_size and query_gap < ref_gap:
                events.append(_detect_deletion(left, right, min_event_size))
            elif query_gap >= min_event_size and ref_gap < query_gap:
                events.append(_detect_insertion(left, right, min_event_size))

    return events


def rearrangement_segment_stream(
    segments: Iterator["pysam.AlignedSegment"],
    max_breakpoint_distance: int = 1000,
    min_segment_length: int = 50,
    min_event_size: int = 30,
    reference_lengths: Optional[Dict[str, int]] = None,
    min_embedded_inversion_size: int = 150,
    min_embedded_overlap_fraction: float = 0.5,
    max_junction_slop: int = 20,
) -> Iterator[Tuple[str, List[RearrangementEvent], List["pysam.AlignedSegment"]]]:
    """Group a query_name-sorted pysam stream and infer rearrangements per read."""
    valid_segments = (segment for segment in segments if segment.cigartuples)

    for query_name, group in groupby(valid_segments, key=lambda s: s.query_name):
        segment_list = list(group)
        mappings = [
            (
                segment.reference_name,
                segment.reference_start,
                segment.is_reverse,
                segment.cigartuples,
            )
            for segment in segment_list
        ]
        events = infer_rearrangements(
            mappings,
            max_breakpoint_distance=max_breakpoint_distance,
            min_segment_length=min_segment_length,
            min_event_size=min_event_size,
            reference_lengths=reference_lengths,
            min_embedded_inversion_size=min_embedded_inversion_size,
            min_embedded_overlap_fraction=min_embedded_overlap_fraction,
            max_junction_slop=max_junction_slop,
        )
        yield query_name, events, segment_list
