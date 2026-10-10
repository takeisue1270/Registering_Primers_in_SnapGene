import pytest

from snapgene_primers import dna_format
from tests.helpers import build_dna_bytes


def test_parse_serialize_roundtrip():
    data = build_dna_bytes("acgtacgtaa", circular=True)
    blocks = dna_format.parse_blocks(data)
    assert [b.type_id for b in blocks] == [
        dna_format.HEADER_BLOCK_TYPE,
        dna_format.SEQUENCE_BLOCK_TYPE,
    ]
    assert dna_format.serialize_blocks(blocks) == data


def test_read_sequence_block_circular_flag():
    data = build_dna_bytes("acgtacgtaa", circular=True)
    blocks = dna_format.parse_blocks(data)
    seq_block = dna_format.find_block(blocks, dna_format.SEQUENCE_BLOCK_TYPE)
    info = dna_format.read_sequence_block(seq_block)
    assert info.circular is True
    assert info.sequence == "acgtacgtaa"


def test_read_sequence_block_linear_flag():
    data = build_dna_bytes("acgtacgtaa", circular=False)
    blocks = dna_format.parse_blocks(data)
    seq_block = dna_format.find_block(blocks, dna_format.SEQUENCE_BLOCK_TYPE)
    info = dna_format.read_sequence_block(seq_block)
    assert info.circular is False


def test_rejects_truncated_file():
    data = build_dna_bytes("acgt")
    with pytest.raises(dna_format.DnaFormatError):
        dna_format.parse_blocks(data[:-2])


def test_rejects_missing_header():
    data = build_dna_bytes("acgt")
    # Drop the header block, leaving only the sequence block.
    blocks = dna_format.parse_blocks(data)
    body_only = dna_format.serialize_blocks(blocks[1:])
    with pytest.raises(dna_format.DnaFormatError):
        dna_format.parse_blocks(body_only)


def test_insertion_index_prefers_after_features_block():
    blocks = [
        dna_format.Block(dna_format.HEADER_BLOCK_TYPE, b""),
        dna_format.Block(dna_format.SEQUENCE_BLOCK_TYPE, b"\x01acgt"),
        dna_format.Block(dna_format.ADDITIONAL_PROPERTIES_BLOCK_TYPE, b""),
        dna_format.Block(dna_format.FEATURES_BLOCK_TYPE, b""),
        dna_format.Block(6, b""),  # Notes
    ]
    assert dna_format.insertion_index_for_primers_block(blocks) == 4


def test_insertion_index_falls_back_to_after_sequence():
    blocks = [
        dna_format.Block(dna_format.HEADER_BLOCK_TYPE, b""),
        dna_format.Block(dna_format.SEQUENCE_BLOCK_TYPE, b"\x01acgt"),
        dna_format.Block(6, b""),
    ]
    assert dna_format.insertion_index_for_primers_block(blocks) == 2
