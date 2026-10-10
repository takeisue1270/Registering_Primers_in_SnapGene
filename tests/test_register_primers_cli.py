import openpyxl

import register_primers
from snapgene_primers import dna_format, primers_xml
from tests.helpers import build_dna_bytes

PLASMID_SEQ = "TTTT" + "ACGTACGTAC" + "GGGG"  # primer sits at 0-based 4-13


def _make_workbook(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Oligo"
    ws.append(["Oligo", "Name", "Sequence", "Concat."])
    ws.append(["O01-01", "Fwd", "ACGTACGTAC", None])
    wb.save(path)


def test_file_option_targets_only_the_given_files_and_ignores_the_rest_of_the_dir(tmp_path):
    plasmid_dir = tmp_path / "PlasmidList"
    plasmid_dir.mkdir()
    target = plasmid_dir / "target.dna"
    target.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False))
    untouched = plasmid_dir / "other.dna"
    untouched.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False))

    xlsx = tmp_path / "nucleotides.xlsx"
    _make_workbook(xlsx)

    report = tmp_path / "report.csv"
    exit_code = register_primers.main(
        [
            "--excel",
            str(xlsx),
            "--file",
            str(target),
            "--apply",
            "--report",
            str(report),
        ]
    )
    assert exit_code == 0

    data = target.read_bytes()
    blocks = dna_format.parse_blocks(data)
    pb = dna_format.find_block(blocks, dna_format.PRIMERS_BLOCK_TYPE)
    _, primers = primers_xml.parse_existing_primers(pb.payload.decode("utf-8"))
    assert [p.name for p in primers] == ["Fwd (O01-01)"]

    # The other file in the same directory must be completely untouched.
    assert untouched.read_bytes() == build_dna_bytes(PLASMID_SEQ, circular=False)


def test_missing_file_argument_reports_error(tmp_path, capsys):
    xlsx = tmp_path / "nucleotides.xlsx"
    _make_workbook(xlsx)

    exit_code = register_primers.main(
        ["--excel", str(xlsx), "--file", str(tmp_path / "does_not_exist.dna")]
    )
    assert exit_code == 1
    assert "見つかりません" in capsys.readouterr().err
