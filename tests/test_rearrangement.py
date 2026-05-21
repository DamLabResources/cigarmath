"""Test rearrangement inference from split-read mappings."""

from collections import Counter, namedtuple

import pytest

import cigarmath as cm
from cigarmath.defn import cigarstr2tup
from cigarmath.inference import inferred_query_sequence_length
from cigarmath.rearrangement import (
    format_read_rearrangement_summary,
    infer_rearrangements,
    rearrangement_segment_stream,
    reference_lengths_from_pysam_header,
    RearrangementEvent,
)

HXB2F_LEN = 9086

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
    return inferred_query_sequence_length(cigartuples)


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


def test_detect_ref_wrap_genomic():
    """Late-genome segment then proximal terminus (HIV linear-ref pattern)."""
    read_len = 6700
    mappings = [
        _reverse_split("HXB2F", 2102, query_start=47, aligned_length=6480, read_len=read_len),
        _reverse_split("HXB2F", 54, query_start=6571, aligned_length=126, read_len=read_len),
    ]
    _assert_uniform_read_length(mappings)
    ref_lengths = {"HXB2F": HXB2F_LEN}
    events = infer_rearrangements(
        mappings,
        max_breakpoint_distance=1000,
        min_segment_length=50,
        min_event_size=30,
        reference_lengths=ref_lengths,
    )
    assert len(events) == 1
    assert events[0].event_type == "REF_WRAP"
    assert events[0].reference_size == 54 - 8582


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
    with pysam.AlignmentFile(test_sam, "r") as handle:
        ref_lengths = reference_lengths_from_pysam_header(handle.header)
    stream = sorted(
        cm.io.segment_stream_pysam(test_sam, mode="r"),
        key=lambda segment: segment.query_name,
    )
    results = list(
        rearrangement_segment_stream(stream, reference_lengths=ref_lengths)
    )

    query_names = [query_name for query_name, _, _ in results]
    assert query_names == sorted(query_names)
    assert len(results) == 193

    reads_with_events = [(q, events, segments) for q, events, segments in results if events]
    assert len(reads_with_events) == 49

    event_counts = Counter(event.event_type for _, events, _ in results for event in events)
    assert event_counts == Counter({"INV": 8, "DEL": 38, "REF_WRAP": 5})

    query_name, events, segments = reads_with_events[0]
    event = events[0]
    assert isinstance(segments[event.segment_indices[0]], pysam.AlignedSegment)
    assert segments[event.segment_indices[0]].query_name == query_name
    assert segments[event.segment_indices[1]].query_name == query_name



def test_format_summary_inversion_two_segment():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        _reverse_split(CHR1, 1380, query_start=500, aligned_length=400, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    summary = format_read_rearrangement_summary(
        "inv2", mappings, events, min_segment_length=1
    )
    assert summary == (
        "read inv2 len=1000 segments=2 events=1\n"
        "[Clip] q[0, 100)\n"
        "[S0 +] q[100,500) ref:chr1:1000-1400\n"
        "[S1 -] q[500,900) ref:chr1:1380-1780 | INV ref -20\n"
        "[Clip] q[900, 1000)"
    )


def test_format_summary_deletion():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=500, read_len=read_len),
        _forward_split(CHR1, 11500, query_start=500, aligned_length=500, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    summary = format_read_rearrangement_summary(
        "del", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    assert summary == (
        "read del len=1000 segments=2 events=1\n"
        "[S0 +] q[0,500) ref:chr1:1000-1500\n"
        "[DEL] ref +10000 bp\n"
        "[S1 +] q[500,1000) ref:chr1:11500-12000"
    )


def test_format_summary_insertion():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=500, read_len=read_len),
        _forward_split(CHR1, 1500, query_start=600, aligned_length=400, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    summary = format_read_rearrangement_summary(
        "ins", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    assert summary == (
        "read ins len=1000 segments=2 events=1\n"
        "[S0 +] q[0,500) ref:chr1:1000-1500\n"
        "[INS] query +100 bp\n"
        "[S1 +] q[600,1000) ref:chr1:1500-1900"
    )


def test_format_summary_duplication():
    read_len = 3000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=1000, read_len=read_len),
        _forward_split(CHR1, 1800, query_start=1000, aligned_length=1000, read_len=read_len),
    ]
    events = infer_rearrangements(
        mappings,
        max_breakpoint_distance=100,
        min_segment_length=1,
        min_event_size=30,
    )
    summary = format_read_rearrangement_summary(
        "dup", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    assert summary == (
        "read dup len=3000 segments=2 events=1\n"
        "[S0 +] q[0,1000) ref:chr1:1000-2000\n"
        "[S1 +] q[1000,2000) ref:chr1:1800-2800 | DUP ref -200"
    )


def test_format_summary_translocation():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=400, read_len=read_len),
        _forward_split(CHR2, 80000, query_start=500, aligned_length=400, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    summary = format_read_rearrangement_summary(
        "tra", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    assert summary == (
        "read tra len=1000 segments=2 events=1\n"
        "[S0 +] q[0,400) ref:chr1:1000-1400\n"
        "[S1 +] q[500,900) ref:chr7:80000-80400 | TRA query +100"
    )


def test_format_summary_mixed_events():
    read_len = 2000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=300, read_len=read_len),
        _reverse_split(CHR1, 1280, query_start=300, aligned_length=300, read_len=read_len),
        _forward_split(CHR1, 1560, query_start=600, aligned_length=300, read_len=read_len),
        _forward_split(CHR1, 12000, query_start=900, aligned_length=500, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=30)
    summary = format_read_rearrangement_summary(
        "mixed", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    assert summary == (
        "read mixed len=2000 segments=4 events=3\n"
        "[S0 +] q[0,300) ref:chr1:1000-1300\n"
        "[S1 -] q[300,600) ref:chr1:1280-1580 | INV ref -20\n"
        "[S2 +] q[600,900) ref:chr1:1560-1860 | INV ref -20\n"
        "[DEL] ref +10140 bp\n"
        "[S3 +] q[900,1400) ref:chr1:12000-12500"
    )


def test_format_summary_clip_threshold():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=100, aligned_length=400, read_len=read_len),
        _forward_split(CHR1, 1400, query_start=500, aligned_length=400, read_len=read_len),
    ]
    summary = format_read_rearrangement_summary(
        "clean", mappings, [], min_segment_length=1, show_clip_threshold=30
    )
    assert "[Clip] q[0, 100)" in summary
    assert "[Clip] q[900, 1000)" in summary

    summary_high = format_read_rearrangement_summary(
        "clean", mappings, [], min_segment_length=1, show_clip_threshold=200
    )
    assert "[Clip]" not in summary_high


def test_format_summary_query_order_and_original_indices():
    read_len = 1000
    mappings = [
        _forward_split(CHR1, 1000, query_start=0, aligned_length=300, read_len=read_len),
        _reverse_split(CHR1, 1280, query_start=300, aligned_length=300, read_len=read_len),
        _forward_split(CHR1, 1560, query_start=600, aligned_length=300, read_len=read_len),
    ]
    events = infer_rearrangements(mappings, min_segment_length=1, min_event_size=1)
    summary = format_read_rearrangement_summary(
        "order", mappings, events, min_segment_length=1, show_clip_threshold=1000
    )
    lines = summary.splitlines()
    segment_lines = [line for line in lines if line.startswith("[S")]
    assert segment_lines == [
        "[S0 +] q[0,300) ref:chr1:1000-1300",
        "[S1 -] q[300,600) ref:chr1:1280-1580 | INV ref -20",
        "[S2 +] q[600,900) ref:chr1:1560-1860 | INV ref -20",
    ]

    shuffled = [mappings[2], mappings[0], mappings[1]]
    shuffled_events = infer_rearrangements(shuffled, min_segment_length=1, min_event_size=1)
    shuffled_summary = format_read_rearrangement_summary(
        "shuffled", shuffled, shuffled_events, min_segment_length=1, show_clip_threshold=1000
    )
    shuffled_segment_lines = [
        line for line in shuffled_summary.splitlines() if line.startswith("[S")
    ]
    assert shuffled_segment_lines == [
        "[S1 +] q[0,300) ref:chr1:1000-1300",
        "[S2 -] q[300,600) ref:chr1:1280-1580 | INV ref -20",
        "[S0 +] q[600,900) ref:chr1:1560-1860 | INV ref -20",
    ]


SAM_SUMMARY_ARTIFACT = "tests/test_data/.rearrangement_summaries_test_sam.txt"


def test_write_test_sam_rearrangement_summaries_artifact():
    """Write all read summaries for manual inspection (git-ignored)."""
    pysam = pytest.importorskip("pysam")

    test_sam = "tests/test_data/test.sam"
    with pysam.AlignmentFile(test_sam, "r") as handle:
        ref_lengths = reference_lengths_from_pysam_header(handle.header)
    stream = sorted(
        cm.io.segment_stream_pysam(test_sam, mode="r"),
        key=lambda segment: segment.query_name,
    )

    blocks = []
    for query_name, events, segments in rearrangement_segment_stream(
        stream, reference_lengths=ref_lengths
    ):
        mappings = [
            (
                segment.reference_name,
                segment.reference_start,
                segment.is_reverse,
                segment.cigartuples,
            )
            for segment in segments
        ]
        blocks.append(
            format_read_rearrangement_summary(query_name, mappings, events)
        )

    artifact_path = SAM_SUMMARY_ARTIFACT
    content = "\n\n".join(blocks) + "\n"
    with open(artifact_path, "w", encoding="utf-8") as handle:
        handle.write(content)

    assert content.strip()
    assert content.count("read ") == 193
