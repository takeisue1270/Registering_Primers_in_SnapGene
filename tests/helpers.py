"""Helpers to build small synthetic .dna byte strings for tests.

None of these fixtures come from real lab data -- they're built from scratch
purely to exercise the block format, so this repository never needs to ship
any proprietary sequence data.
"""

from __future__ import annotations

from snapgene_primers import dna_format


def build_dna_bytes(
    sequence: str, *, circular: bool = True, primers_xml_text: str | None = None
) -> bytes:
    blocks = [
        dna_format.Block(dna_format.HEADER_BLOCK_TYPE, b"SnapGene fixture"),
        dna_format.Block(
            dna_format.SEQUENCE_BLOCK_TYPE,
            bytes([dna_format.TOPOLOGY_CIRCULAR_BIT if circular else 0])
            + sequence.encode("ascii"),
        ),
    ]
    if primers_xml_text is not None:
        blocks.append(dna_format.Block(dna_format.PRIMERS_BLOCK_TYPE, primers_xml_text.encode("utf-8")))
    return dna_format.serialize_blocks(blocks)
