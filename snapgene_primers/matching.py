"""Exact-match search of an oligo sequence against a (possibly circular) plasmid."""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _find_all(haystack: str, needle: str) -> list[int]:
    """0-based start indices of every (overlapping-allowed) exact match, case-insensitive."""
    if not needle:
        return []
    hay = haystack.lower()
    ndl = needle.lower()
    positions = []
    start = 0
    while True:
        idx = hay.find(ndl, start)
        if idx == -1:
            return positions
        positions.append(idx)
        start = idx + 1


@dataclass(frozen=True)
class MatchLocation:
    """A single exact binding site.

    ``start0``/``end0`` are 0-based, inclusive, in SnapGene's own Primer
    ``BindingSite location="start-end"`` coordinate system (verified against
    713 pre-existing exact-match primer entries in real .dna files; see
    tests/test_matching.py). For a site that wraps past the origin of a
    circular plasmid, ``end0 < start0`` -- this mirrors what SnapGene itself
    writes for such sites.

    ``bound_strand`` is ``1`` when the primer anneals to the plasmid's top
    (plus) strand -- i.e. the primer's own sequence equals the
    reverse-complement of the plasmid at this location -- and ``0`` when the
    primer anneals to the bottom (minus) strand, i.e. the primer's sequence
    equals the plasmid's top-strand bases directly.
    """

    start0: int
    end0: int
    bound_strand: int

    @property
    def location(self) -> str:
        return f"{self.start0}-{self.end0}"

    @property
    def wraps_origin(self) -> bool:
        return self.end0 < self.start0


def expected_bases(plasmid_seq: str, match: MatchLocation, primer_len: int) -> str:
    """Re-extract the plasmid bases spanned by ``match`` (for verification).

    Returns the bases exactly as they appear in ``plasmid_seq`` (whatever
    case SnapGene stored), reading through the circular origin if needed.
    """
    n = len(plasmid_seq)
    if not match.wraps_origin:
        return plasmid_seq[match.start0 : match.start0 + primer_len]
    return plasmid_seq[match.start0 :] + plasmid_seq[: match.end0 + 1]


def find_primer_matches(
    plasmid_seq: str, circular: bool, primer_seq: str
) -> list[MatchLocation]:
    """Find every exact, full-length binding site of ``primer_seq`` in ``plasmid_seq``.

    Both orientations are searched (primer literally present on the top
    strand, and primer present as the reverse complement of the top strand).
    For circular plasmids, sites that wrap past the origin are also found.
    """
    n = len(plasmid_seq)
    plen = len(primer_seq)
    if plen == 0 or plen > n:
        return []

    search_space = plasmid_seq
    if circular and plen > 1:
        search_space = plasmid_seq + plasmid_seq[: plen - 1]

    seen: set[tuple[int, int, int]] = set()
    matches: list[MatchLocation] = []

    def add(idx0: int, bound_strand: int) -> None:
        if idx0 >= n:
            # Only occurrences fully described starting within the original
            # sequence are kept; matches starting in the duplicated tail are
            # the same physical site as one already found starting earlier.
            return
        start0 = idx0
        end0 = (idx0 + plen - 1) % n
        key = (start0, end0, bound_strand)
        if key in seen:
            return
        seen.add(key)
        matches.append(MatchLocation(start0, end0, bound_strand))

    for idx in _find_all(search_space, primer_seq):
        add(idx, bound_strand=0)

    rc = reverse_complement(primer_seq)
    for idx in _find_all(search_space, rc):
        add(idx, bound_strand=1)

    return matches
