from snapgene_primers import dna_format
from snapgene_primers.oligo_source import OligoRecord
from snapgene_primers.primers_xml import parse_existing_primers
from snapgene_primers.registration import apply_file, plan_file

from tests.helpers import build_dna_bytes

PLASMID_SEQ = "TTTT" + "ACGTACGTAC" + "GGGG" + "CCCCAAAACC" + "TTTT"
#               0-3     4-13           14-17   18-27          28-31


def test_plan_finds_forward_match_and_creates_primer_when_no_block_exists(tmp_path):
    dna_path = tmp_path / "plasmid.dna"
    dna_path.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=None))

    oligo = OligoRecord(oligo_id="O01-01", name="Fwd", sequence="ACGTACGTAC", row=2)
    result = plan_file(str(dna_path), [oligo])

    assert result.error is None
    assert result.already_registered == []
    assert len(result.to_register) == 1
    entry = result.to_register[0]
    assert entry.matches[0].location == "4-13"
    assert entry.matches[0].bound_strand == 0

    apply_file(result, backup_dir=str(tmp_path / "backup"))
    assert result.applied is True
    assert result.backup_path is not None

    data = dna_path.read_bytes()
    blocks = dna_format.parse_blocks(data)
    primers_block = dna_format.find_block(blocks, dna_format.PRIMERS_BLOCK_TYPE)
    assert primers_block is not None
    _, primers = parse_existing_primers(primers_block.payload.decode("utf-8"))
    assert len(primers) == 1
    assert primers[0].name == "Fwd (O01-01)"
    assert primers[0].sequence == "ACGTACGTAC"


def test_rerun_is_idempotent_and_skips_already_registered(tmp_path):
    dna_path = tmp_path / "plasmid.dna"
    dna_path.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=None))
    oligo = OligoRecord(oligo_id="O01-01", name="Fwd", sequence="ACGTACGTAC", row=2)

    first = plan_file(str(dna_path), [oligo])
    apply_file(first, backup_dir=str(tmp_path / "backup"))

    second = plan_file(str(dna_path), [oligo])
    assert second.to_register == []
    assert len(second.already_registered) == 1
    assert second.already_registered[0].oligo_id == "O01-01"


def test_no_match_oligo_is_neither_registered_nor_flagged_as_duplicate(tmp_path):
    dna_path = tmp_path / "plasmid.dna"
    dna_path.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=None))
    oligo = OligoRecord(oligo_id="O99-99", name="NoMatch", sequence="GGGGGGGGGGGG", row=2)

    result = plan_file(str(dna_path), [oligo])
    assert result.to_register == []
    assert result.already_registered == []


def test_reverse_complement_match_is_registered_with_bound_strand_one(tmp_path):
    dna_path = tmp_path / "plasmid.dna"
    dna_path.write_bytes(build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=None))
    # Reverse complement of the 4-13 segment "ACGTACGTAC" is "GTACGTACGT".
    oligo = OligoRecord(oligo_id="O01-02", name="Rev", sequence="GTACGTACGT", row=3)

    result = plan_file(str(dna_path), [oligo])
    assert len(result.to_register) == 1
    match = result.to_register[0].matches[0]
    assert match.location == "4-13"
    assert match.bound_strand == 1

    apply_file(result, backup_dir=str(tmp_path / "backup"))
    data = dna_path.read_bytes()
    blocks = dna_format.parse_blocks(data)
    primers_block = dna_format.find_block(blocks, dna_format.PRIMERS_BLOCK_TYPE)
    xml_text = primers_block.payload.decode("utf-8")
    assert 'boundStrand="1"' in xml_text
    assert 'location="4-13"' in xml_text


def test_circular_wraparound_primer_is_registered(tmp_path):
    seq = "CCTTGG"
    dna_path = tmp_path / "small_circular.dna"
    dna_path.write_bytes(build_dna_bytes(seq, circular=True, primers_xml_text=None))
    oligo = OligoRecord(oligo_id="O05-01", name="Wrap", sequence="TGGCC", row=2)

    result = plan_file(str(dna_path), [oligo])
    assert len(result.to_register) == 1
    assert result.to_register[0].matches[0].location == "3-1"

    apply_file(result, backup_dir=str(tmp_path / "backup"))
    data = dna_path.read_bytes()
    blocks = dna_format.parse_blocks(data)
    primers_block = dna_format.find_block(blocks, dna_format.PRIMERS_BLOCK_TYPE)
    assert 'location="3-1"' in primers_block.payload.decode("utf-8")


def test_apply_backs_up_original_file_before_writing(tmp_path):
    dna_path = tmp_path / "plasmid.dna"
    original_bytes = build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=None)
    dna_path.write_bytes(original_bytes)
    oligo = OligoRecord(oligo_id="O01-01", name="Fwd", sequence="ACGTACGTAC", row=2)

    result = plan_file(str(dna_path), [oligo])
    backup_dir = tmp_path / "backup"
    apply_file(result, backup_dir=str(backup_dir))

    backup_file = backup_dir / "plasmid.dna"
    assert backup_file.exists()
    assert backup_file.read_bytes() == original_bytes
    assert dna_path.read_bytes() != original_bytes


def test_existing_primer_with_same_oligo_id_is_detected_as_duplicate(tmp_path):
    existing_xml = (
        '<?xml version="1.0"?><Primers nextValidID="1">'
        '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
        'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
        'minimumFivePrimeAnnealing="15"/>'
        '<Primer recentID="0" name="Fwd (O01-01)" sequence="ACGTACGTAC" '
        'description="&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;" '
        'dateAdded="2026-01-01T00:00:00Z">'
        '<BindingSite location="4-13" boundStrand="0" annealedBases="ACGTACGTAC" '
        'meltingTemperature="55">'
        '<Component hybridizedRange="4-13" bases="ACGTACGTAC"/></BindingSite>'
        "</Primer></Primers>"
    )
    dna_path = tmp_path / "plasmid.dna"
    dna_path.write_bytes(
        build_dna_bytes(PLASMID_SEQ, circular=False, primers_xml_text=existing_xml)
    )
    # Same oligo ID, even with a different display name -> still detected as already registered.
    oligo = OligoRecord(oligo_id="O01-01", name="RenamedFwd", sequence="ACGTACGTAC", row=2)

    result = plan_file(str(dna_path), [oligo])
    assert result.to_register == []
    assert len(result.already_registered) == 1
