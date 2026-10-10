# snapgene-primer-registrar

*[日本語版 README はこちら](README.ja.md)*

Matches oligo sequences from an Excel oligo-management sheet against SnapGene
`.dna` plasmid files and automatically registers, as a "primer" on the
plasmid, only the oligos whose sequence is a **full, exact match**.

It's meant to be re-run every time new oligos are added to the sheet:
primers that are already registered are detected and skipped automatically
(idempotent).

## What it does

- Scans every `.dna` file in a given folder
- For each oligo sequence listed in the `Oligo` sheet of the workbook, searches for a **full, exact match** on either strand, including binding sites that wrap past the origin of a circular plasmid
- Only on a match, writes a `<Primer>` entry into the `.dna` file in the same format SnapGene itself uses when you add a primer through the GUI
- The primer is named `Name (ID)`, rebuilt from those two columns — the same string the sheet's `Concat.` column shows, so the name matches the sheet whether or not that column is filled in
- Skips an oligo in a given file if a primer containing that same oligo ID (e.g. `O01-01`) is already registered there
- Defaults to a **dry run** that writes nothing and only produces a report; pass `--apply` to actually write
- Automatically backs up each file before writing, and re-reads what actually landed on disk afterwards to independently re-verify position and sequence (see "Safety design" below)

## Install

```powershell
pip install -r requirements.txt
```

## The sheet it reads

The workbook is only ever read, never written. It has to be laid out like this,
because the tool addresses one column positionally and the others by name:

| Where | Must contain | Example |
| --- | --- | --- |
| Row 1 | The header row | |
| Column A | The oligo ID. Addressed **by position**, because this column is headed `Oligo` as a section label rather than a column name. | `O01-01` |
| A column headed `Name` | The oligo name | `Example_Del.F` |
| A column headed `Sequence` | The sequence, A/C/G/T only, case-insensitive | `ACGTACGT…` |

Anything else in the sheet is ignored, and the `Name` and `Sequence` columns may
sit in any order at any position. A row is skipped silently when its ID, name or
sequence is empty, or when the sequence contains anything other than A, C, G or
T — degenerate bases and modifications included. The run stops with an error if
the sheet named by `--sheet` does not exist, or if it has no `Name` or
`Sequence` header.

## Usage

Make sure the workbook is **not currently open in Excel** (it will fail to read while file-locked).

```powershell
# Dry run first: nothing is written, only a report is produced
python register_primers.py --excel "path\to\Nucleotide.List.xlsx" --plasmid-dir "path\to\PlasmidList"

# Review the report, then actually write with --apply
python register_primers.py --excel "path\to\Nucleotide.List.xlsx" --plasmid-dir "path\to\PlasmidList" --apply
```

`--excel` and `--plasmid-dir` default to `./InSilico/Nucleotide.List.xlsx` and
`./InSilico/Vector/PlasmidList`, so if you run from a project root laid out
that way you can omit them entirely:

```powershell
python register_primers.py            # dry run
python register_primers.py --apply    # actually register
```

Whenever you add new oligos, just re-run the second line above. Anything
already registered is skipped automatically; only the new ones get added.

To register into one or two specific `.dna` files instead of scanning a whole
folder (for example, reference plasmids that live outside `--plasmid-dir`),
pass `--file` once per file and skip `--plasmid-dir` entirely:

```powershell
python register_primers.py --file "path\to\LocusA.dna" --file "path\to\LocusB.dna" --apply
```

### Main options

| Option | Description | Default |
| --- | --- | --- |
| `--excel` | Nucleotide List workbook (reads the `Oligo` sheet) | `./InSilico/Nucleotide.List.xlsx` |
| `--sheet` | Sheet name | `Oligo` |
| `--plasmid-dir` | Folder of `.dna` files, searched recursively | `./InSilico/Vector/PlasmidList` |
| `--file` | Register into this specific `.dna` file instead of scanning `--plasmid-dir`. Repeatable. | (none) |
| `--apply` | Actually write changes (omit for a dry run) | off |
| `--backup-dir` | Where to back up files before writing | a sibling of `<plasmid-dir>` (or of the first `--file`'s folder) named `_primer_registration_backups/<timestamp>` |
| `--report` | CSV report path | `./primer_registration_report_<timestamp>.csv` |

## Safety design

Because this tool writes directly into `.dna` files that may be in active use
in the lab, it follows these rules:

1. **Dry run by default**: nothing is changed unless `--apply` is given.
2. **Automatic backup**: before writing, the pre-write bytes of each modified file are copied as-is to `--backup-dir`.
3. **Atomic write**: writes go to a temp file first, then `os.replace()` swaps it in, so a crash mid-write can never corrupt the original.
4. **Double position verification**, the part that matters most:
   - Before writing: the bases at the computed binding-site location (0-based, against the plasmid sequence) are re-extracted and checked against the oligo sequence being registered (or its reverse complement, if it binds the opposite strand). The write only proceeds if they match exactly; otherwise an exception aborts that entry.
   - After writing: the file is re-read from disk, and the newly added primer's name, sequence, and binding-site position are independently re-verified against what is actually on disk. Any mismatch triggers an immediate restore from the just-made backup.
   - This design mechanically guarantees that a computed base position always matches the actual base at that position in the real plasmid sequence.
5. **CSV report**: every run writes out which oligo was registered in which file at which position, and which were already registered and skipped, so results stay auditable over time.

## About the `.dna` file format

SnapGene's `.dna` format is a binary container made of blocks, each a 1-byte
type + a 4-byte (big-endian) length + payload. This tool reverse-engineered
and verified its coordinate system against 713 primer entries that SnapGene
itself had written into existing `.dna` files. The cases that pinned the
coordinate convention down are reproduced as synthetic fixtures in
`tests/test_matching.py`.

- A primer's binding site location is `BindingSite location="start-end"`, **0-based and inclusive on both ends**.
- `boundStrand="0"` means the primer sequence appears literally on the plasmid's plus strand (i.e. it anneals to the minus strand).
- `boundStrand="1"` means the primer sequence is the reverse complement of the plasmid's plus strand at that location (i.e. it anneals to the plus strand).
- A binding site that wraps past the origin of a circular plasmid is written with `end < start` (e.g. `"2808-49"`).

New primers are added by treating the existing `Primers` block's XML as text:
a new `<Primer>...</Primer>` fragment is appended and `nextValidID` is
updated, without touching the bytes of any existing entry. This "text
surgery" approach minimizes the risk of corrupting existing data (see
`snapgene_primers/primers_xml.py` for details).

## Testing

```powershell
pip install -r requirements-dev.txt
pytest
```

All tests use only synthetic data built in the test code; no experimental data is included.

## Limitations

- Only full, exact matches (100% identity, case-insensitive) are registered; partial matches that tolerate mismatches are not supported.
- If the same oligo matches multiple locations within the same plasmid, all sites found are registered as multiple `BindingSite` entries under a single `<Primer>` element.
- The displayed Tm (melting temperature) is a simple approximation (Wallace rule / GC% formula) and does not exactly match SnapGene's internal calculation; it is never used in any matching/verification logic.
- The `Oligo` sheet in the workbook is only ever read, never modified.

## The format is not a documented one

`.dna` is a proprietary format with no published specification. Everything this
tool knows about it was worked out by reading files SnapGene had written, as
described above. Two consequences follow, and you should accept both before
pointing it at files you care about:

- **A SnapGene update can break it.** If the format changes, this tool may
  write a file that SnapGene then reads incorrectly or refuses to open. The
  backup taken before every write is your way back; keep it until you have
  opened the file in SnapGene and checked it.
- **Nothing here is endorsed by anyone.** SnapGene is a commercial product of
  Dotmatics, and "SnapGene" is their trademark. This project is independent:
  not affiliated with, endorsed by, or supported by them, and it uses the name
  only to say which files it reads.

Open a few written files in SnapGene the first time you use this on a new
version, and read the CSV report rather than assuming the run went well.

## Disclaimer

This is a research tool, provided as is and without warranty of any kind, as
the [licence](LICENSE) states. It is not intended for diagnostic, therapeutic
or clinical use.

It modifies files you depend on. The dry run, the backup, and the verification
before and after each write exist because of that, but they do not transfer the
responsibility: the annotations it produces are yours to check, and a primer
registered at a verified position is still only as correct as the sequence in
the file and the row in the sheet. Nothing here substitutes for confirming a
primer experimentally.

## Citation

[CITATION.cff](CITATION.cff) holds the metadata, so GitHub's **Cite this
repository** button works.

