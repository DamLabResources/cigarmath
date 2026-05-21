"""Infer chromosomal rearrangements from split-read alignments."""

__copyright__ = """Copyright (C) 2022-present
    Dampier & DV Klopfenstein, PhD.
    All rights reserved"""
__author__ = "Will Dampier, PhD"

from collections import namedtuple
from dataclasses import dataclass
from itertools import groupby
from typing import Iterator, List, Optional, Tuple, TYPE_CHECKING

from cigarmath.block import query_block, reference_block
from cigarmath.clipping import left_clipping, right_clipping
from cigarmath.defn import CigarTuples

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


def _normalize_segments(
    mappings: List[Mapping],
    min_segment_length: int,
) -> List[_Segment]:
    if not mappings:
        return []

    read_len = max(
        query_block(cigartuples)[1] + right_clipping(cigartuples)
        for _, _, _, cigartuples in mappings
    )

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
    return segments


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


def _detect_duplication(
    left: _Segment,
    right: _Segment,
    min_event_size: int,
) -> Optional[RearrangementEvent]:
    reference_size = right.ref_start - left.ref_end
    if reference_size >= 0:
        return None
    if abs(reference_size) < min_event_size:
        return None
    query_gap = max(0, right.q_start - left.q_end)
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
) -> List[RearrangementEvent]:
    """Classify rearrangements between adjacent split-read mapping segments."""
    segments = _normalize_segments(mappings, min_segment_length)
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
            event = _detect_duplication(left, right, min_event_size)
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
        )
        yield query_name, events, segment_list
