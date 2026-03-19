"""Domain and interaction site extractor.

Extracts binding sites, active sites, and interaction regions from
UniProt feature annotations, slices the canonical sequence to produce
fragment sequences, and tags each with its ligand/partner description.
"""

from __future__ import annotations

from core.models import SequenceFeature


BINDING_TYPES = {"binding site"}
ACTIVE_TYPES = {"active site"}
DOMAIN_TYPES = {"domain", "region", "motif", "zinc finger", "coiled coil", "transmembrane"}

INTERACTION_KEYWORDS = frozenset(
    ["interaction", "binding", "receptor", "ligand", "interface", "contact"]
)


def _make_feature(feat_dict: dict) -> SequenceFeature:
    return SequenceFeature(**feat_dict)


def extract_features_from_uniprot(
    features: list[dict], canonical_sequence: str
) -> tuple[list[SequenceFeature], list[SequenceFeature], list[SequenceFeature]]:
    """Parse UniProt feature list and return (binding_sites, active_sites, domains).

    Args:
        features:           Raw list of feature dicts from UniProt JSON response.
        canonical_sequence: Full amino acid sequence string.

    Returns:
        Tuple of three lists: binding_sites, active_sites, domains.
    """
    binding_sites: list[SequenceFeature] = []
    active_sites: list[SequenceFeature] = []
    domains: list[SequenceFeature] = []

    for feat in features:
        ftype_raw = feat.get("type", "")
        ftype = ftype_raw.lower()
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")

        if start is None or end is None:
            continue

        fragment = canonical_sequence[start - 1: end] if canonical_sequence else ""
        desc = feat.get("description") or ""
        ligand_raw = feat.get("ligand")
        if isinstance(ligand_raw, dict):
            ligand = ligand_raw.get("name", "") or ""
        elif isinstance(ligand_raw, str):
            ligand = ligand_raw
        else:
            ligand = ""

        sf = SequenceFeature(
            feature_type=ftype_raw,
            start=start,
            end=end,
            sequence_fragment=fragment,
            description=desc or None,
            ligand=ligand or None,
        )

        if ftype in BINDING_TYPES:
            binding_sites.append(sf)
        elif ftype in ACTIVE_TYPES:
            active_sites.append(sf)
        elif ftype in DOMAIN_TYPES:
            # Also catch regions with interaction keywords
            if ftype == "region":
                if any(kw in desc.lower() for kw in INTERACTION_KEYWORDS):
                    binding_sites.append(sf)
                else:
                    domains.append(sf)
            else:
                domains.append(sf)

    return binding_sites, active_sites, domains


def format_feature_summary(sf: SequenceFeature) -> str:
    """Format a SequenceFeature for human-readable display."""
    lines = [
        f"Feature type:  {sf.feature_type}",
        f"Position:      {sf.start}-{sf.end}",
        f"Sequence:      {sf.sequence_fragment}",
    ]
    if sf.ligand:
        lines.append(f"Ligand:        {sf.ligand}")
    if sf.description:
        lines.append(f"Description:   {sf.description}")
    return "\n".join(lines)
