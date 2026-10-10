import openpyxl
import pytest

from snapgene_primers.oligo_source import OligoSheetError, load_oligos


def _make_workbook(path, rows):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Oligo"
    ws.append(["Oligo", "Name", "Sequence", "Concat."])
    for row in rows:
        ws.append(row)
    wb.save(path)


def test_loads_valid_rows_and_skips_incomplete_ones(tmp_path):
    xlsx = tmp_path / "nucleotides.xlsx"
    _make_workbook(
        xlsx,
        [
            ["O01-01", "Example_Del.F", "ACGTACGTACGTTTAAGCATTAGCATAGCA", None],
            ["O01-02", "Example_173.R", "ttaaggccttaaggccttaaggcc", None],
            ["O01-03", "NoSequence", None, None],  # missing sequence -> skipped
            [None, "NoId", "ACGT", None],  # missing id -> skipped
        ],
    )
    records = load_oligos(str(xlsx))
    assert [r.oligo_id for r in records] == ["O01-01", "O01-02"]
    assert records[0].name == "Example_Del.F"
    assert records[0].sequence == "ACGTACGTACGTTTAAGCATTAGCATAGCA"
    assert records[0].concat_name == "Example_Del.F (O01-01)"


def test_skips_sequences_with_invalid_bases(tmp_path):
    xlsx = tmp_path / "nucleotides.xlsx"
    _make_workbook(xlsx, [["O01-01", "Bad", "ACGXT", None]])
    assert load_oligos(str(xlsx)) == []


def test_missing_sheet_raises(tmp_path):
    xlsx = tmp_path / "nucleotides.xlsx"
    _make_workbook(xlsx, [["O01-01", "Ok", "ACGT", None]])
    with pytest.raises(OligoSheetError):
        load_oligos(str(xlsx), sheet_name="DoesNotExist")


def test_missing_required_column_raises(tmp_path):
    xlsx = tmp_path / "nucleotides.xlsx"
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Oligo"
    ws.append(["Oligo", "NotSequence"])
    ws.append(["O01-01", "x"])
    wb.save(xlsx)
    with pytest.raises(OligoSheetError):
        load_oligos(str(xlsx))
