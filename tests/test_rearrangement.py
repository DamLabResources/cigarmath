"""Test rearrangement inference from split-read mappings."""

from collections import Counter, namedtuple

import pytest

import cigarmath as cm
from cigarmath.block import query_block
from cigarmath.clipping import right_clipping
from cigarmath.defn import cigarstr2tup
from cigarmath.rearrangement import (
    infer_rearrangements,
    rearrangement_segment_stream,
    RearrangementEvent,
)

CHR1 = "chr1"
CHR2 = "chr7"

FakeSegment = namedtuple(
    "FakeSegment",
    [
        "query_name",
        "reference_name",
        "reference_start",
        "is_reverse",
        "cigartuples",
    ],
)


def _inferred_query_length(cigartuples):
    """Total read length implied by one supplementary alignment record."""
    return query_block(cigartuples)[1] + right_clipping(cigartuples)


def _mapping(ref_name, ref_start, is_reverse, cigar):
    return (ref_name, ref_start, is_reverse, cigarstr2tup(cigar))


def _assert_uniform_read_length(mappings):
    lengths = [_inferred_query_length(cigartuples) for *_, cigartuples in mappings]
    assert len(set(lengths)) == 1, f"inconsistent read lengths: {lengths}"


def _forward_split(ref_name, ref_start, query_start, aligned_length, read_len):
    """Build a forward-strand split mapping covering query [query_start, query_end)."""
    query_end = query_start + aligned_length
    left = query_start
    right = read_len - query_end
    cigar = f"{left}S {aligned_length}M {right}S"
    return _mapping(ref_name, ref_start, False, cigar)


def _reverse_split(ref_name, ref_start, query_start, aligned_length, read_len):
    """Build a reverse-strand split mapping projecting to the same query interval."""
    query_end = query_start + aligned_length
    left = read_len - query_end
    right = query_start
    cigar = f"{left}S {aligned_length}M {right}S"
    return _mapping(ref_name, ref_start, True, cigar)


def test_no_event_single_segment():
    # ref_start=1000, 0S400M600S
    mapping = _forward_split(CHR1, 1000, query_start=0, aligned_length=400, read_len=1000)
    assert _inferred_query_length(mapping[3]) == 1000
    events = infer_rearrangements([mapping])
    assert events == []


def test_no_event_clean_split():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _forward_split(CHR1, 1400, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=1400, 500S 400M 100S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    assert events == []


def test_detect_inversion_two_segment():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _reverse_split(CHR1, 1380, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=1380, 100S 400M 500S (reverse)
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    assert len(events) == 1
    assert events[0].event_type == "INV"
    assert events[0].segment_indices == (0, 1)


def test_detect_inversion_three_segment():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=300, read_len=read_len),
        # ref_start=1000, 0S 300M 700S
        _reverse_split(CHR1, 1280, query_start=300, aligned_length=300, read_len=read_len),
        # ref_start=1280, 700S 300M 300S (reverse)
        _forward_split(CHR1, 1560, query_start=600, aligned_length=300, read_len=read_len),
        # ref_start=1560, 600S 300M 100S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    assert len(events) == 2
    assert all(event.event_type == "INV" for event in events)


def test_detect_deletion():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=500, read_len=read_len),
        # ref_start=1000, 0S 500M 500S
        _forward_split(CHR1, 11500, query_start=500, aligned_length=500, read_len=read_len),
        # ref_start=11500, 500S 500M 0S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    assert len(events) == 1
    assert events[0].event_type == "DEL"
    assert events[0].reference_size == 10000
    assert events[0].query_size == 0


def test_detect_insertion():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=500, read_len=read_len),
        # ref_start=1000, 0S 500M 500S
        _forward_split(CHR1, 1500, query_start=600, aligned_length=400, read_len=read_len),
        # ref_start=1500, 600S 400M 0S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    assert len(events) == 1
    assert events[0].event_type == "INS"
    assert events[0].query_size == 100
    assert events[0].reference_size == 0


def test_detect_duplication():
    read_len = 3000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=1000, read_len=read_len),
        # ref_start=1000, 0S 1000M 2000S
        _forward_split(CHR1, 1800, query_start=1000, aligned_length=1000, read_len=read_len),
        # ref_start=1800, 1000S 1000M 1000S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(
        mappings,
        max_breakpoint_distance=100,
        min_segment_length=1,
        min_event_size=30,
    )
    assert len(events) == 1
    assert events[0].event_type == "DUP"
    assert events[0].reference_size == -200


def test_detect_translocation():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=400, read_len=read_len),
        # ref_start=1000, 0S 400M 600S
        _forward_split(CHR2, 80000, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=80000, 500S 400M 100S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    assert len(events) == 1
    assert events[0].event_type == "TRA"
    assert events[0].reference_size is None


def test_input_order_independence():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=300, read_len=read_len),
        # ref_start=1000, 0S 300M 700S
        _reverse_split(CHR1, 1280, query_start=300, aligned_length=300, read_len=read_len),
        # ref_start=1280, 700S 300M 300S (reverse)
        _forward_split(CHR1, 1560, query_start=600, aligned_length=300, read_len=read_len),
        # ref_start=1560, 600S 300M 100S
    ]
    _assert_uniform_read_length(mappings)
    expected = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    shuffled = [mappings[2], mappings[0], mappings[1]]
    actual = infer_rearrangements(shuffled, min_segment_length=1, min_event_size=1)
    assert len(actual) == len(expected)
    assert {event.event_type for event in actual} == {event.event_type for event in expected}
    assert len(actual) == 2


def test_min_segment_length_filters_noise():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _reverse_split(CHR1, 1390, query_start=500, aligned_length=10, read_len=read_len),
        # ref_start=1390, 490S 10M 500S (reverse)
        _forward_split(CHR1, 50000, query_start=510, aligned_length=400, read_len=read_len),
        # ref_start=50000, 510S 400M 90S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(
        mappings,
        min_segment_length=50,
        min_event_size=1,
        max_breakpoint_distance=1000,
    )
    assert all(event.event_type != "INV" for event in events)


def test_min_event_size_filters_small_indels():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _forward_split(CHR1, 1405, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=1405, 500S 400M 100S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    assert events == []


def test_mixed_events():
    read_len = 2000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=300, read_len=read_len),
        # ref_start=1000, 0S 300M 1700S
        _reverse_split(CHR1, 1280, query_start=300, aligned_length=300, read_len=read_len),
        # ref_start=1280, 1700S 300M 300S (reverse)
        _forward_split(CHR1, 1560, query_start=600, aligned_length=300, read_len=read_len),
        # ref_start=1560, 600S 300M 1100S
        _forward_split(CHR1, 12000, query_start=900, aligned_length=500, read_len=read_len),
        # ref_start=12000, 900S 500M 600S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    assert {event.event_type for event in events} == {"INV", "DEL"}
    assert sum(event.event_type == "INV" for event in events) == 2


def test_event_repr_is_human_readable():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=500, read_len=read_len),
        # ref_start=1000, 0S 500M 500S
        _forward_split(CHR1, 11500, query_start=500, aligned_length=500, read_len=read_len),
        # ref_start=11500, 500S 500M 0S
    ]
    _assert_uniform_read_length(mappings)
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    text = str(events[0])
    assert "DEL" in text
    assert CHR1 in text
    assert "1000" in text or "1500" in text


def _fake_segment(query_name, mapping):
    ref_name, ref_start, is_reverse, cigartuples = mapping
    return FakeSegment(query_name, ref_name, ref_start, is_reverse, cigartuples)


def test_rearrangement_segment_stream_groups_by_query_name():
    read_len = 1000
    read_a = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _reverse_split(CHR1, 1380, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=1380, 100S 400M 500S (reverse)
    ]
    read_b = [
        _forward_split(CHR1, 2000, query_start=0, aligned_length=400, read_len=read_len),
        # ref_start=2000, 0S400M 600S
    ]
    _assert_uniform_read_length(read_a)
    _assert_uniform_read_length(read_b)
    segments = [
        _fake_segment("read_a", read_a[0]),
        _fake_segment("read_a", read_a[1]),
        _fake_segment("read_b", read_b[0]),
    ]
    results = list(
        rearrangement_segment_stream(iter(segments), min_segment_length=1, min_event_size=1)
    )
    assert len(results) == 2
    assert results[0][0] == "read_a"
    assert len(results[0][1]) == 1
    assert results[1][0] == "read_b"
    assert results[1][1] == []


def test_rearrangement_segment_stream_empty_events_for_single_segment():
    # ref_start=1000, 0S400M600S
    mapping = _forward_split(CHR1, 1000, query_start=0, aligned_length=400, read_len=1000)
    segment = _fake_segment("solo", mapping)
    results = list(
        rearrangement_segment_stream(
            iter([segment]),
            min_segment_length=1,
            min_event_size=1,
        )
    )
    assert results == [("solo", [], [segment])]


def test_rearrangement_segment_stream_emits_events():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        # ref_start=1000, 100S 400M 500S
        _reverse_split(CHR1, 1380, query_start=500, aligned_length=400, read_len=read_len),
        # ref_start=1380, 100S 400M 500S (reverse)
    ]
    _assert_uniform_read_length(mappings)
    segments = [_fake_segment("inv_read", mapping) for mapping in mappings]
    query_name, events, yielded = next(
        rearrangement_segment_stream(iter(segments), min_segment_length=1, min_event_size=1)
    )
    assert query_name == "inv_read"
    assert len(events) == 1
    assert events[0].segment_indices == (0, 1)
    assert yielded[events[0].segment_indices[0]] is segments[0]
    assert yielded[events[0].segment_indices[1]] is segments[1]


def test_rearrangement_segment_stream_skips_missing_cigar():
    # ref_start=1380, 100S 400M 500S (reverse)
    mapping = _reverse_split(CHR1, 1380, query_start=500, aligned_length=400, read_len=1000)
    segments = [
        FakeSegment("read_x", CHR1, 1000, False, None),
        _fake_segment("read_x", mapping),
    ]
    results = list(
        rearrangement_segment_stream(iter(segments), min_segment_length=1, min_event_size=1)
    )
    assert len(results) == 1
    assert results[0][0] == "read_x"
    assert results[0][1] == []


def test_rearrangement_segment_stream_test_sam():
    """End-to-end on tests/test_data/test.sam (query_name-sorted)."""
    pysam = pytest.importorskip("pysam")

    test_sam = "tests/test_data/test.sam"
    stream = sorted(
        cm.io.segment_stream_pysam(test_sam, mode="r"),
        key=lambda segment: segment.query_name,
    )
    results = list(rearrangement_segment_stream(stream))

    query_names = [query_name for query_name, _, _ in results]
    assert query_names == sorted(query_names)
    assert len(results) == 193

    reads_with_events = [(q, events, segments) for q, events, segments in results if events]
    assert len(reads_with_events) == 49

    event_counts = Counter(event.event_type for _, events, _ in results for event in events)
    assert event_counts == Counter({"INV": 8, "DEL": 38, "DUP": 5})

    for query_name, events, segments in reads_with_events:
        print(query_name, events, segments)

    query_name, events, segments = reads_with_events[0]
    event = events[0]
    assert isinstance(segments[event.segment_indices[0]], pysam.AlignedSegment)
    assert segments[event.segment_indices[0]].query_name == query_name
    assert segments[event.segment_indices[1]].query_name == query_name

    for _, _, segments in results:
        if len(segments) < 2:
            continue
        lengths = [
            query_block(segment.cigartuples)[1] + right_clipping(segment.cigartuples)
            for segment in segments
            if segment.cigartuples
        ]
        assert len(set(lengths)) == 1, f"{segments[0].query_name}: {lengths}"
