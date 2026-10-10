"""Low-level read/write support for the SnapGene ``.dna`` binary container format.

A ``.dna`` file is a flat sequence of length-prefixed blocks::

    +----------+------------------+-----------------------------+
    | type (1B)| payload length   | payload (length bytes)      |
    |          | (4B, big-endian) |                              |
    +----------+------------------+-----------------------------+

Block type ``0`` holds the DNA sequence itself (1 flag byte + ASCII bases).
Block type ``5`` holds the "Primers" table as an XML fragment.  Everything
else is left untouched by this tool.
"""

from __future__ import annotations

from dataclasses import dataclass

HEADER_BLOCK_TYPE = 9
SEQUENCE_BLOCK_TYPE = 0
ADDITIONAL_PROPERTIES_BLOCK_TYPE = 8
FEATURES_BLOCK_TYPE = 10
PRIMERS_BLOCK_TYPE = 5

# Bit 0 of the sequence-block flag byte marks the topology as circular.
TOPOLOGY_CIRCULAR_BIT = 0x01


@dataclass
class Block:
    type_id: int
    payload: bytes


class DnaFormatError(ValueError):
    """Raised when a ``.dna`` file does not look like a valid SnapGene container."""


def parse_blocks(data: bytes) -> list[Block]:
    blocks: list[Block] = []
    i, n = 0, len(data)
    while i < n:
        if i + 5 > n:
            raise DnaFormatError(f"truncated block header at offset {i}")
        type_id = data[i]
        length = int.from_bytes(data[i + 1 : i + 5], "big")
        start = i + 5
        end = start + length
        if end > n:
            raise DnaFormatError(f"truncated block payload at offset {i} (type {type_id})")
        blocks.append(Block(type_id, data[start:end]))
        i = end
    if not blocks or blocks[0].type_id != HEADER_BLOCK_TYPE:
        raise DnaFormatError("file does not start with a SnapGene cookie block (type 9)")
    return blocks


def serialize_blocks(blocks: list[Block]) -> bytes:
    out = bytearray()
    for block in blocks:
        out.append(block.type_id)
        out += len(block.payload).to_bytes(4, "big")
        out += block.payload
    return bytes(out)


@dataclass
class SequenceInfo:
    circular: bool
    sequence: str  # exactly as stored in the file (case as SnapGene wrote it)


def read_sequence_block(block: Block) -> SequenceInfo:
    if not block.payload:
        raise DnaFormatError("sequence block is empty")
    flags = block.payload[0]
    sequence = block.payload[1:].decode("ascii")
    return SequenceInfo(circular=bool(flags & TOPOLOGY_CIRCULAR_BIT), sequence=sequence)


def find_block(blocks: list[Block], type_id: int) -> Block | None:
    for block in blocks:
        if block.type_id == type_id:
            return block
    return None


def insertion_index_for_primers_block(blocks: list[Block]) -> int:
    """Where to insert a brand-new Primers block if one does not already exist.

    Observed SnapGene ordering is ``... Additional Properties(8), Features(10),
    Primers(5), Notes(6) ...``. We insert right after the last Features/
    Additional-Properties block, falling back to right after the sequence
    block, so a freshly created block lands in a sane place even though
    SnapGene does not actually require any particular block order.
    """
    last_anchor = None
    for idx, block in enumerate(blocks):
        if block.type_id in (ADDITIONAL_PROPERTIES_BLOCK_TYPE, FEATURES_BLOCK_TYPE):
            last_anchor = idx
    if last_anchor is not None:
        return last_anchor + 1
    for idx, block in enumerate(blocks):
        if block.type_id == SEQUENCE_BLOCK_TYPE:
            return idx + 1
    return len(blocks)
