from snapgene_primers.matching import (
    expected_bases,
    find_primer_matches,
    reverse_complement,
)


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACCGGTT") == "AACCGGTT"
    assert reverse_complement("acgtACGT") == "ACGTacgt"


def test_forward_strand_match_is_bound_strand_zero():
    plasmid = "TTTT" + "ACGTACGTAC" + "TTTT"
    matches = find_primer_matches(plasmid, circular=False, primer_seq="ACGTACGTAC")
    assert matches == [_ml(4, 13, 0)]


def test_reverse_complement_match_is_bound_strand_one():
    plasmid = "TTTT" + "ACGTACGTAC" + "TTTT"
    # This primer's reverse complement is "ACGTACGTAC", the segment embedded above.
    primer = "GTACGTACGT"
    assert reverse_complement(primer) == "ACGTACGTAC"
    matches = find_primer_matches(plasmid, circular=False, primer_seq=primer)
    assert matches == [_ml(4, 13, 1)]


def test_case_insensitive_match():
    plasmid = "tttt" + "ACGTACGTAC" + "tttt"
    matches = find_primer_matches(plasmid, circular=False, primer_seq="acgtacgtac")
    assert matches == [_ml(4, 13, 0)]


def test_no_match_returns_empty_list():
    plasmid = "AAAAAAAAAAAAAAAA"
    assert find_primer_matches(plasmid, circular=False, primer_seq="GGGGGGGG") == []


def test_primer_longer_than_plasmid():
    assert find_primer_matches("ACGT", circular=True, primer_seq="ACGTACGTACGT") == []


def test_circular_wraparound_match_forward():
    plasmid = "CCTTGG"  # 0:C 1:C 2:T 3:T 4:G 5:G, circular
    matches = find_primer_matches(plasmid, circular=True, primer_seq="TGGCC")
    assert matches == [_ml(3, 1, 0)]
    assert matches[0].wraps_origin is True


def test_wraparound_not_found_when_linear():
    plasmid = "CCTTGG"
    matches = find_primer_matches(plasmid, circular=False, primer_seq="TGGCC")
    assert matches == []


def test_expected_bases_non_wrapping():
    plasmid = "TTTTACGTACGTACTTTT"
    match = _ml(4, 13, 0)
    assert expected_bases(plasmid, match, 10) == "ACGTACGTAC"


def test_expected_bases_wrapping_matches_plasmid_bases():
    plasmid = "CCTTGG"
    match = _ml(3, 1, 0)
    assert expected_bases(plasmid, match, 5) == "TGGCC"


def test_palindrome_matches_both_strands_at_same_location():
    plasmid = "TTTT" + "GAATTC" + "TTTT"  # EcoRI site, palindromic
    matches = find_primer_matches(plasmid, circular=False, primer_seq="GAATTC")
    locations = {(m.start0, m.end0, m.bound_strand) for m in matches}
    assert locations == {(4, 9, 0), (4, 9, 1)}


def _ml(start0, end0, bound_strand):
    from snapgene_primers.matching import MatchLocation

    return MatchLocation(start0, end0, bound_strand)
