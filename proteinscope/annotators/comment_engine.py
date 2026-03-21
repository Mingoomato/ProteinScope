"""Comment engine — rule-based comment generation from ProteinRecord data.

All functions are pure/synchronous and produce deterministic text from
data already present in the ProteinRecord. No AI calls here.
"""

from __future__ import annotations

from typing import Optional

from core.models import ProteinRecord


MAX_ARROW_CHARS = 42


def build_interaction_comment(node: dict, record: ProteinRecord) -> str:
    """Build a sentence describing the protein's role in this pathway node.

    Priority: metabolic reaction data > signaling pathway data > fallback.
    """
    gene = node["gene"]

    # Check metabolic reactions
    for mp in record.metabolic_pathways:
        for rxn in mp.reactions:
            if rxn.equation and gene.upper() in (record.gene_name.upper(),):
                subs = ", ".join(rxn.substrates[:2]) or "substrate"
                prods = ", ".join(rxn.products[:2]) or "product"
                cofactors = ", ".join(rxn.cofactors_required[:2])
                cof_str = f" (cofactor: {cofactors})" if cofactors else ""
                return (f"{gene} catalyzes conversion of {subs} to {prods}{cof_str}")

    # Check signaling pathways
    for sp in record.signaling_pathways:
        if gene.upper() in sp.pathway_name.upper() or gene.upper() == record.gene_name.upper():
            acts = ", ".join(sp.activates[:3]) or "—"
            inh  = ", ".join(sp.inhibits[:3]) or "—"
            role = sp.protein_role or "signaling component"
            return (f"{gene} acts as {role}; activates: {acts}; inhibits: {inh}")

    return f"{gene} participates in this pathway"


def build_drug_comment(node: dict, record: ProteinRecord) -> Optional[str]:
    """Build a drug interaction comment if drug data exists for this gene."""
    gene = node["gene"]
    # Match by gene symbol OR by the queried protein (most common case)
    relevant = [
        d for d in record.drug_interactions
        if gene.upper() == record.gene_name.upper() or
           gene.upper() in [s.upper() for s in d.sources]
    ]
    if not relevant:
        return None
    lines = []
    for d in relevant[:3]:
        ev = d.clinical_significance or "evidence unspecified"
        lines.append(f"{d.drug_name} ({d.interaction_type}, {ev})")
    return "Drug targets: " + "; ".join(lines)


def build_disease_comment(node: dict, record: ProteinRecord) -> Optional[str]:
    """Build a disease association comment if disease data exists."""
    gene = node["gene"]
    names = []
    for dp in record.disease_pathways:
        if (gene.upper() in dp.disease_name.upper() or
                gene.upper() == record.gene_name.upper() or
                any(gene.upper() in s.event.upper() for s in dp.cascade)):
            names.append(dp.disease_name)
    if not names:
        return None
    unique_names = list(dict.fromkeys(names))[:2]
    return "Disease association: " + ", ".join(unique_names)


def build_stream_comment(
    node: dict,
    record: ProteinRecord,
) -> tuple[Optional[str], Optional[str]]:
    """Build upstream and downstream comments from signaling pathway data.

    Returns (upstream_comment, downstream_comment).
    """
    gene = node["gene"]
    upstream = downstream = None

    for sp in record.signaling_pathways:
        if gene.upper() != record.gene_name.upper():
            break
        if sp.activated_by:
            upstream = "Activated by: " + ", ".join(sp.activated_by[:3])
        if sp.activates:
            downstream = "Downstream targets: " + ", ".join(sp.activates[:3])
        if upstream or downstream:
            break

    return upstream, downstream


def build_arrow_label(interaction_comment: str,
                       drug_comment: Optional[str]) -> str:
    """Build a short inline label for the arrow (max MAX_ARROW_CHARS chars).

    Prioritizes drug information (most clinically salient).
    """
    if drug_comment:
        # Extract first drug name only
        first_drug = drug_comment.replace("Drug targets: ", "").split(";")[0].strip()
        label = f"Drug target: {first_drug}"
    else:
        # Use first ~7 words of interaction comment
        words = interaction_comment.split()
        label = " ".join(words[:7])

    if len(label) > MAX_ARROW_CHARS:
        label = label[:MAX_ARROW_CHARS - 3] + "..."
    return label


def get_protein_name(gene_symbol: str, record: ProteinRecord) -> str:
    """Return the full protein name for a gene symbol (use record's name if same gene)."""
    if gene_symbol.upper() == record.gene_name.upper():
        return record.protein_name
    return gene_symbol  # fallback: just return the symbol
