import xml.etree.ElementTree as ET

from snapgene_primers.primers_xml import (
    BindingSiteSpec,
    empty_primers_block,
    insert_primers,
    parse_existing_primers,
    render_primer_element,
)

_SAMPLE_BLOCK = (
    '<?xml version="1.0"?><Primers nextValidID="2">'
    '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
    'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
    'minimumFivePrimeAnnealing="15"/>'
    '<Primer recentID="0" name="ExistingFwd (O01-01)" sequence="ACGTACGTAC" '
    'description="&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;" '
    'dateAdded="2026-01-01T00:00:00Z">'
    '<BindingSite location="4-13" boundStrand="0" annealedBases="ACGTACGTAC" '
    'meltingTemperature="55">'
    '<Component hybridizedRange="4-13" bases="ACGTACGTAC"/></BindingSite>'
    '<BindingSite simplified="1" location="4-13" boundStrand="0" '
    'annealedBases="ACGTACGTAC" meltingTemperature="55">'
    '<Component hybridizedRange="4-13" bases="ACGTACGTAC"/></BindingSite>'
    "</Primer></Primers>"
)


def test_parse_existing_primers():
    next_id, primers = parse_existing_primers(_SAMPLE_BLOCK)
    assert next_id == 2
    assert len(primers) == 1
    assert primers[0].recent_id == 0
    assert primers[0].name == "ExistingFwd (O01-01)"
    assert primers[0].sequence == "ACGTACGTAC"


def test_render_primer_element_is_well_formed_and_roundtrips():
    site = BindingSiteSpec(
        location="4-13", bound_strand=0, annealed_bases="ACGTACGTAC", melting_temperature=55
    )
    fragment = render_primer_element(
        recent_id=1,
        name='New "Primer" & Friend (O02-02)',
        sequence="ACGTACGTAC",
        date_added="2026-02-02T00:00:00Z",
        binding_sites=[site],
    )
    wrapped = f"<Primers>{fragment}</Primers>"
    root = ET.fromstring(wrapped)  # raises if not well-formed
    primer_el = root.find("Primer")
    assert primer_el.get("name") == 'New "Primer" & Friend (O02-02)'
    assert primer_el.get("sequence") == "ACGTACGTAC"
    binding_sites = primer_el.findall("BindingSite")
    assert len(binding_sites) == 2  # one plain + one simplified
    assert binding_sites[0].get("location") == "4-13"
    assert binding_sites[1].get("simplified") == "1"


def test_insert_primers_appends_and_bumps_next_valid_id():
    fragment = render_primer_element(
        recent_id=2,
        name="NewOne (O02-02)",
        sequence="TTTTGGGG",
        date_added="2026-02-02T00:00:00Z",
        binding_sites=[
            BindingSiteSpec(
                location="0-7", bound_strand=0, annealed_bases="TTTTGGGG", melting_temperature=40
            )
        ],
    )
    updated = insert_primers(_SAMPLE_BLOCK, [fragment], new_next_valid_id=3)

    assert 'nextValidID="3"' in updated
    assert "ExistingFwd (O01-01)" in updated  # original entry untouched
    assert "NewOne (O02-02)" in updated

    next_id, primers = parse_existing_primers(updated)
    assert next_id == 3
    names = {p.name for p in primers}
    assert names == {"ExistingFwd (O01-01)", "NewOne (O02-02)"}


def test_insert_into_empty_primers_block():
    fragment = render_primer_element(
        recent_id=0,
        name="First (O99-01)",
        sequence="AAAACCCC",
        date_added="2026-03-03T00:00:00Z",
        binding_sites=[
            BindingSiteSpec(
                location="0-7", bound_strand=0, annealed_bases="AAAACCCC", melting_temperature=38
            )
        ],
    )
    updated = insert_primers(empty_primers_block(), [fragment], new_next_valid_id=1)
    next_id, primers = parse_existing_primers(updated)
    assert next_id == 1
    assert primers[0].name == "First (O99-01)"
