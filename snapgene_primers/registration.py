"""Plan and (optionally) apply primer registration into a directory of .dna files."""

from __future__ import annotations

import datetime as dt
import os
import re
import shutil
from dataclasses import dataclass, field

from . import dna_format, matching, primers_xml
from .melting_temp import estimate_tm
from .oligo_source import OligoRecord


class VerificationError(RuntimeError):
    """Raised when a computed binding-site position does not reproduce the
    expected sequence when re-extracted from the plasmid. This should never
    trigger if matching.py is correct; it exists as a hard safety net so a
    logic bug can never silently write a wrong primer position."""


@dataclass
class RegistrationEntry:
    oligo: OligoRecord
    matches: list[matching.MatchLocation]


@dataclass
class FileResult:
    dna_path: str
    circular: bool = False
    seq_len: int = 0
    to_register: list[RegistrationEntry] = field(default_factory=list)
    already_registered: list[OligoRecord] = field(default_factory=list)
    error: str | None = None
    applied: bool = False
    backup_path: str | None = None


def _oligo_id_token(oligo_id: str) -> re.Pattern[str]:
    return re.compile(r"\(" + re.escape(oligo_id) + r"\)")


def _load_sequence_and_existing(
    data: bytes,
) -> tuple[dna_format.SequenceInfo, list, list["dna_format.Block"]]:
    blocks = dna_format.parse_blocks(data)
    seq_block = dna_format.find_block(blocks, dna_format.SEQUENCE_BLOCK_TYPE)
    if seq_block is None:
        raise VerificationError("DNA配列ブロック(type 0)が見つかりません")
    seq_info = dna_format.read_sequence_block(seq_block)

    primers_block = dna_format.find_block(blocks, dna_format.PRIMERS_BLOCK_TYPE)
    if primers_block is not None:
        _, existing = primers_xml.parse_existing_primers(
            primers_block.payload.decode("utf-8")
        )
    else:
        existing = []
    return seq_info, existing, blocks


def plan_file(dna_path: str, oligos: list[OligoRecord]) -> FileResult:
    """Read-only pass: decide, for each oligo, whether it is already
    registered in this file, matches and should be newly registered, or
    doesn't match at all (in which case it is simply omitted)."""
    try:
        data = open(dna_path, "rb").read()
        seq_info, existing, _blocks = _load_sequence_and_existing(data)
    except (dna_format.DnaFormatError, VerificationError) as exc:
        return FileResult(dna_path, error=str(exc))

    result = FileResult(dna_path, circular=seq_info.circular, seq_len=len(seq_info.sequence))

    for oligo in oligos:
        token = _oligo_id_token(oligo.oligo_id)
        if any(token.search(p.name) for p in existing):
            result.already_registered.append(oligo)
            continue
        matches = matching.find_primer_matches(
            seq_info.sequence, seq_info.circular, oligo.sequence
        )
        if matches:
            _verify_matches(seq_info.sequence, oligo, matches)
            result.to_register.append(RegistrationEntry(oligo, matches))

    return result


def _verify_matches(
    plasmid_seq: str, oligo: OligoRecord, matches: list[matching.MatchLocation]
) -> None:
    """Re-extract the bases at each computed location and confirm they
    really do equal the oligo's sequence (or its reverse complement,
    depending on strand) before we ever trust that location."""
    plen = len(oligo.sequence)
    for match in matches:
        bases = matching.expected_bases(plasmid_seq, match, plen)
        expected = (
            matching.reverse_complement(oligo.sequence)
            if match.bound_strand == 1
            else oligo.sequence
        )
        if bases.lower() != expected.lower():
            raise VerificationError(
                f"位置検証エラー: {oligo.concat_name} の {match.location} "
                f"(boundStrand={match.bound_strand}) を再抽出した配列 '{bases}' が "
                f"期待値 '{expected}' と一致しません。登録を中止します。"
            )


def apply_file(result: FileResult, *, backup_dir: str) -> None:
    """Write ``result.to_register`` into the .dna file on disk.

    Safety sequence: re-read + re-verify from disk -> back up the original
    -> atomic write -> re-read the file that actually landed on disk and
    verify it *again*, restoring the backup if anything is wrong.
    """
    if not result.to_register:
        return

    data = open(result.dna_path, "rb").read()
    seq_info, _existing, blocks = _load_sequence_and_existing(data)

    for entry in result.to_register:
        entry.matches = matching.find_primer_matches(
            seq_info.sequence, seq_info.circular, entry.oligo.sequence
        )
        if not entry.matches:
            raise VerificationError(
                f"{entry.oligo.concat_name}: 書込み直前の再検索でマッチが消失しました "
                "(ファイルが計画時から変更された可能性があります)"
            )
        _verify_matches(seq_info.sequence, entry.oligo, entry.matches)

    primers_idx = next(
        (i for i, b in enumerate(blocks) if b.type_id == dna_format.PRIMERS_BLOCK_TYPE),
        None,
    )
    xml_text = (
        blocks[primers_idx].payload.decode("utf-8")
        if primers_idx is not None
        else primers_xml.empty_primers_block()
    )
    next_valid_id, _ = primers_xml.parse_existing_primers(xml_text)

    date_added = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    fragments = []
    for entry in result.to_register:
        sites = [
            primers_xml.BindingSiteSpec(
                location=m.location,
                bound_strand=m.bound_strand,
                annealed_bases=entry.oligo.sequence,
                melting_temperature=estimate_tm(entry.oligo.sequence),
            )
            for m in entry.matches
        ]
        fragments.append(
            primers_xml.render_primer_element(
                recent_id=next_valid_id,
                name=entry.oligo.concat_name,
                sequence=entry.oligo.sequence,
                date_added=date_added,
                binding_sites=sites,
            )
        )
        next_valid_id += 1

    new_xml = primers_xml.insert_primers(xml_text, fragments, next_valid_id)
    new_block = dna_format.Block(dna_format.PRIMERS_BLOCK_TYPE, new_xml.encode("utf-8"))

    if primers_idx is not None:
        blocks[primers_idx] = new_block
    else:
        blocks.insert(dna_format.insertion_index_for_primers_block(blocks), new_block)

    new_data = dna_format.serialize_blocks(blocks)

    os.makedirs(backup_dir, exist_ok=True)
    backup_path = os.path.join(backup_dir, os.path.basename(result.dna_path))
    shutil.copy2(result.dna_path, backup_path)
    result.backup_path = backup_path

    tmp_path = result.dna_path + ".tmp_write"
    try:
        with open(tmp_path, "wb") as fh:
            fh.write(new_data)
        os.replace(tmp_path, result.dna_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    try:
        _verify_written_file(result)
    except Exception:
        shutil.copy2(backup_path, result.dna_path)
        raise

    result.applied = True


def _verify_written_file(result: FileResult) -> None:
    """Independently re-read the bytes that actually landed on disk (not
    anything held in memory) and confirm every newly written primer's name,
    sequence, and binding-site position round-trip correctly."""
    data = open(result.dna_path, "rb").read()
    seq_info, existing, _blocks = _load_sequence_and_existing(data)
    by_name = {p.name: p for p in existing}

    for entry in result.to_register:
        primer = by_name.get(entry.oligo.concat_name)
        if primer is None:
            raise VerificationError(
                f"書込み後の検証失敗: '{entry.oligo.concat_name}' が読み戻せません"
            )
        if primer.sequence != entry.oligo.sequence:
            raise VerificationError(
                f"書込み後の検証失敗: '{entry.oligo.concat_name}' の配列が一致しません "
                f"(書込み後='{primer.sequence}' / 期待='{entry.oligo.sequence}')"
            )
        _verify_matches(seq_info.sequence, entry.oligo, entry.matches)
