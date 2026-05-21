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
from cigarmath.defn import CigarTuples
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


def _segment_line(segment: _Segment, event_suffix: str = "") -> str:
    line = (
        f"[S{segment.index} {_strand_symbol(segment.is_reverse)}] "
        f"q[{segment.q_start},{segment.q_end}) "
        f"ref:{segment.ref_name}:{segment.ref_start}-{segment.ref_end}"
    )
    if event_suffix:
        line = f"{line} {event_suffix}"
    return line


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
    downstream_labels: dict[int, List[str]] = {}
    for event in events:
        if event.event_type in ("INS", "DEL"):
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
        lines.append(_segment_line(segment, suffix))

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
    return min(
        abs(right.ref_start - left.ref_end),
        abs(left.ref_start - right.ref_end),
    )


def _make_event(
    event_type: str,
    left: _Segment,
    right: _Segment,
    reference_size: Optional[int],
    query_size: int,
) -> RearrangementEvent:
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
    query_gap = max(0, right.q_start - left.q_end)
    return _make_event("TRA", left, right, None, query_gap)


def _detect_inversion(
    left: _Segment,
    right: _Segment,
    max_breakpoint_distance: int,
) -> Optional[RearrangementEvent]:
    if _reference_proximity(left, right) > max_breakpoint_distance:
        return None
    query_gap = max(0, right.q_start - left.q_end)
    return _make_event("INV", left, right, right.ref_start - left.ref_end, query_gap)


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
    """True when a negative ref jump is genome-scale (linear ref wrap), not tandem DUP."""
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
) -> List[RearrangementEvent]:
    """Classify rearrangements between adjacent split-read mapping segments."""
    _, segments = _normalize_segments(mappings, min_segment_length)
    if len(segments) < 2:
        return []

    events: List[RearrangementEvent] = []
    for left, right in zip(segments, segments[1:]):
        if left.ref_name != right.ref_name:
            events.append(_detect_translocation(left, right))
        elif left.is_reverse != right.is_reverse:
            event = _detect_inversion(left, right, max_breakpoint_distance)
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
        )
        yield query_name, events, segment_list
