#!/usr/bin/env python3
"""Register oligos as SnapGene primers wherever their sequence matches exactly.

Scans every ``.dna`` file under --plasmid-dir; for each oligo listed in the
"Oligo" sheet of --excel, registers it as a primer in any plasmid where its
sequence appears as a full, exact match (either strand, including across the
origin of a circular plasmid). Oligos already registered in a given file
(matched by their oligo ID, e.g. "O01-01") are left untouched, so this script
is safe -- and intended -- to re-run every time new oligos are added.

By default this is a dry run: nothing is written, only a report is produced.
Pass --apply to actually modify files (each modified file is backed up first,
and every write is verified by re-reading it back from disk).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import sys
from pathlib import Path

from snapgene_primers.oligo_source import OligoSheetError, load_oligos
from snapgene_primers.registration import FileResult, apply_file, plan_file


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--excel",
        default="./InSilico/Nucleotide.List.xlsx",
        help="Nucleotide List workbook path (default: %(default)s)",
    )
    parser.add_argument(
        "--sheet", default="Oligo", help="Sheet name to read oligos from (default: %(default)s)"
    )
    parser.add_argument(
        "--plasmid-dir",
        default="./InSilico/Vector/PlasmidList",
        help="Directory of .dna plasmid files, searched recursively (default: %(default)s)",
    )
    parser.add_argument(
        "--file",
        dest="files",
        action="append",
        metavar="PATH",
        help=(
            "Register primers into this specific .dna file instead of scanning "
            "--plasmid-dir. Repeatable to target more than one file."
        ),
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually write changes. Without this flag, only a report is produced.",
    )
    parser.add_argument(
        "--backup-dir",
        default=None,
        help=(
            "Where to back up each modified .dna file before writing. Must not be "
            "inside --plasmid-dir, or its own backup copies would be picked up as "
            "input on the next run "
            "(default: a sibling of <plasmid-dir> named _primer_registration_backups/<timestamp>)"
        ),
    )
    parser.add_argument(
        "--report",
        default=None,
        help="CSV report path (default: ./primer_registration_report_<timestamp>.csv)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    try:
        oligos = load_oligos(args.excel, sheet_name=args.sheet)
    except OligoSheetError as exc:
        print(f"エラー: {exc}", file=sys.stderr)
        return 1
    print(f"{len(oligos)} 件のオリゴ配列を '{args.excel}' の '{args.sheet}' シートから読み込みました。")

    if args.files:
        dna_files = [Path(f).resolve() for f in args.files]
        missing = [f for f in dna_files if not f.is_file()]
        if missing:
            print(
                "エラー: 次の --file で指定されたファイルが見つかりません: "
                + ", ".join(str(m) for m in missing),
                file=sys.stderr,
            )
            return 1
        backup_dir = Path(args.backup_dir).resolve() if args.backup_dir else (
            dna_files[0].parent.parent / _BACKUP_DIR_NAME / timestamp
        )
        print(f"{len(dna_files)} 件の指定ファイルを対象にします: " + ", ".join(f.name for f in dna_files))
    else:
        plasmid_dir = Path(args.plasmid_dir).resolve()
        backup_dir = Path(args.backup_dir).resolve() if args.backup_dir else (
            plasmid_dir.parent / _BACKUP_DIR_NAME / timestamp
        )
        if backup_dir == plasmid_dir or backup_dir in plasmid_dir.parents:
            print(
                f"エラー: --backup-dir ('{backup_dir}') が --plasmid-dir ('{plasmid_dir}') "
                "と同じか、その親ディレクトリになっています。指定し直してください。",
                file=sys.stderr,
            )
            return 1

        all_dna_files = sorted(plasmid_dir.rglob("*.dna"))
        dna_files = discover_plasmid_files(all_dna_files, backup_dir)
        skipped_backup_files = len(all_dna_files) - len(dna_files)
        if not dna_files:
            print(f"エラー: '{plasmid_dir}' に .dna ファイルが見つかりません。", file=sys.stderr)
            return 1
        print(f"{len(dna_files)} 件の .dna ファイルを '{plasmid_dir}' から検出しました。")
        if skipped_backup_files:
            print(f"  (バックアップフォルダ内の {skipped_backup_files} 件は対象から除外しました)")

    report_path = args.report or f"primer_registration_report_{timestamp}.csv"

    results: list[FileResult] = []
    plan_errors = 0
    for dna_path in dna_files:
        result = plan_file(str(dna_path), oligos)
        if result.error:
            plan_errors += 1
            print(f"  [解析エラー] {dna_path.name}: {result.error}", file=sys.stderr)
        results.append(result)

    matched_oligo_ids = {e.oligo.oligo_id for r in results for e in r.to_register}
    matched_oligo_ids |= {o.oligo_id for r in results for o in r.already_registered}
    unmatched_count = len(oligos) - len(matched_oligo_ids)

    total_new = sum(len(r.to_register) for r in results)
    total_dup = sum(len(r.already_registered) for r in results)
    print(
        f"\n新規登録予定: {total_new} 件 / 登録済みでスキップ: {total_dup} 件 / "
        f"どのプラスミドにも一致しなかったオリゴ: {unmatched_count} 件"
    )

    _write_report(report_path, results)
    print(f"レポートを書き出しました: {report_path}")

    if not args.apply:
        print("\n--apply が指定されていないため、ドライランのみでファイルは変更していません。")
        return 1 if plan_errors else 0

    apply_errors = 0
    for result in results:
        if not result.to_register:
            continue
        try:
            apply_file(result, backup_dir=str(backup_dir))
            names = ", ".join(e.oligo.concat_name for e in result.to_register)
            print(f"  [登録完了] {Path(result.dna_path).name}: {names}")
        except Exception as exc:  # noqa: BLE001 - surfaced to the operator, one file must not abort the rest
            apply_errors += 1
            print(f"  [登録失敗] {Path(result.dna_path).name}: {exc}", file=sys.stderr)

    if any(r.applied for r in results):
        print(f"\n変更前のバックアップ: {backup_dir}")

    return 1 if (plan_errors or apply_errors) else 0


_BACKUP_DIR_NAME = "_primer_registration_backups"


def discover_plasmid_files(all_dna_files: list[Path], backup_dir: Path) -> list[Path]:
    """Filter out anything under ``backup_dir`` -- or under any directory
    literally named ``_primer_registration_backups`` wherever it appears --
    so that a previous run's backups are never mistaken for plasmid input on
    a later run (see tests/test_register_primers.py for the regression this
    guards against)."""
    return [
        p
        for p in all_dna_files
        if not _is_under(p, backup_dir) and _BACKUP_DIR_NAME not in p.parts
    ]


def _is_under(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def _write_report(report_path: str, results: list[FileResult]) -> None:
    with open(report_path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["dna_file", "action", "oligo_id", "oligo_name", "primer_name", "locations"]
        )
        for result in results:
            for entry in result.to_register:
                locations = "; ".join(
                    f"{m.location}(strand={m.bound_strand})" for m in entry.matches
                )
                writer.writerow(
                    [
                        Path(result.dna_path).name,
                        "register",
                        entry.oligo.oligo_id,
                        entry.oligo.name,
                        entry.oligo.concat_name,
                        locations,
                    ]
                )
            for oligo in result.already_registered:
                writer.writerow(
                    [
                        Path(result.dna_path).name,
                        "already_registered",
                        oligo.oligo_id,
                        oligo.name,
                        oligo.concat_name,
                        "",
                    ]
                )


if __name__ == "__main__":
    raise SystemExit(main())
