"""Read/write the ``Primers`` XML fragment stored in a .dna block-5 payload.

New primers are spliced into the existing XML text rather than round-tripped
through a generic XML serializer, so byte-for-byte content that SnapGene (or
a human) already wrote -- attribute order, entity encoding, whitespace -- is
left completely untouched. Only two things are ever changed on an existing
block: the ``nextValidID`` counter, and the insertion of new ``<Primer>``
elements immediately before the closing ``</Primers>`` tag.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from xml.sax.saxutils import quoteattr
import xml.etree.ElementTree as ET

_EMPTY_TEMPLATE = (
    '<?xml version="1.0"?><Primers nextValidID="0">'
    '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
    'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
    'minimumFivePrimeAnnealing="15"/></Primers>'
)

_NEXT_VALID_ID_RE = re.compile(r'(<Primers\b[^>]*\bnextValidID=")(\d+)(")')
_CLOSING_TAG_RE = re.compile(r"</Primers>\s*$")


@dataclass
class ExistingPrimer:
    recent_id: int
    name: str
    sequence: str


class PrimersBlockError(ValueError):
    pass


def parse_existing_primers(xml_text: str) -> tuple[int, list[ExistingPrimer]]:
    """Return (next_valid_id, [existing primers]) from a Primers block's XML text."""
    root = ET.fromstring(xml_text)
    if root.tag != "Primers":
        raise PrimersBlockError(f"unexpected root element <{root.tag}>")
    next_valid_id = int(root.get("nextValidID", "0"))
    primers = []
    for el in root.findall("Primer"):
        recent_id = int(el.get("recentID", "0"))
        primers.append(
            ExistingPrimer(
                recent_id=recent_id,
                name=el.get("name", ""),
                sequence=el.get("sequence", ""),
            )
        )
    return next_valid_id, primers


@dataclass
class BindingSiteSpec:
    location: str
    bound_strand: int
    annealed_bases: str
    melting_temperature: int


def _render_binding_site(spec: BindingSiteSpec, *, simplified: bool) -> str:
    simplified_attr = ' simplified="1"' if simplified else ""
    return (
        f'<BindingSite{simplified_attr} location="{spec.location}" '
        f'boundStrand="{spec.bound_strand}" '
        f"annealedBases={quoteattr(spec.annealed_bases)} "
        f'meltingTemperature="{spec.melting_temperature}">'
        f'<Component hybridizedRange="{spec.location}" '
        f"bases={quoteattr(spec.annealed_bases)}/>"
        "</BindingSite>"
    )


def render_primer_element(
    *,
    recent_id: int,
    name: str,
    sequence: str,
    date_added: str,
    binding_sites: list[BindingSiteSpec],
) -> str:
    sites_xml = "".join(_render_binding_site(s, simplified=False) for s in binding_sites)
    sites_xml += "".join(_render_binding_site(s, simplified=True) for s in binding_sites)
    return (
        f'<Primer recentID="{recent_id}" name={quoteattr(name)} '
        f"sequence={quoteattr(sequence)} "
        'description="&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;" '
        f'dateAdded="{date_added}">{sites_xml}</Primer>'
    )


def insert_primers(xml_text: str, new_primer_fragments: list[str], new_next_valid_id: int) -> str:
    """Return ``xml_text`` with ``new_primer_fragments`` appended and the
    ``nextValidID`` counter updated, without disturbing anything else."""
    if not _CLOSING_TAG_RE.search(xml_text):
        raise PrimersBlockError("could not find closing </Primers> tag to splice into")
    updated = _NEXT_VALID_ID_RE.sub(rf"\g<1>{new_next_valid_id}\g<3>", xml_text, count=1)
    if updated == xml_text:
        raise PrimersBlockError("could not find nextValidID attribute to update")
    insertion = "".join(new_primer_fragments)
    updated = _CLOSING_TAG_RE.sub(insertion + "</Primers>", updated, count=1)
    return updated


def empty_primers_block() -> str:
    return _EMPTY_TEMPLATE
