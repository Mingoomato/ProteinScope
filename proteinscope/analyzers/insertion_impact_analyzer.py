"""Metabolic insertion impact analyzer.

When a foreign gene/protein is inserted into a host organism, this module answers:
  - Which host metabolic pathways would be disrupted?
  - Which host enzymes compete for the same substrates or are inhibited by products?
  - Which host proteins would physically interact with the inserted protein?
  - What cofactors would be depleted and which host pathways depend on them?

Analysis pipeline:
1. Resolve source EC number (UniProt → BRENDA → manual)
2. Fetch KEGG EC entry → catalyzed reaction, substrates, products, cofactors
3. Find human KEGG pathways containing the same EC number → competing enzymes
4. Find human KEGG pathways whose metabolites overlap with source substrates/products
5. STRING-DB: predict host PPI partners (ortholog-transferred interactions)
6. Reactome: fetch source protein's pathway context
7. MetaCyc: cross-check enzyme reactions in source organism
8. Gemini 2.5 Pro: bioinformatics synthesis with conflict prioritisation
"""

from __future__ import annotations

import asyncio
import json as _json
import re
from typing import Optional

import httpx
from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class PathwayConflict(BaseModel):
    pathway_id: str
    pathway_name: str
    conflict_type: str          # substrate_competition | product_inhibition |
                                # cofactor_depletion | ppi_sequestration |
                                # pathway_bypass | product_toxicity
    host_genes_affected: list[str] = Field(default_factory=list)
    description: str
    severity: str               # High | Moderate | Low
    evidence_sources: list[str] = Field(default_factory=list)


class InsertionImpactReport(BaseModel):
    source_gene: str
    source_organism: str
    target_organism: str
    source_ec_number: Optional[str] = None

    # Reaction chemistry of the inserted enzyme
    reaction_equation: Optional[str] = None
    substrates: list[str] = Field(default_factory=list)
    products: list[str] = Field(default_factory=list)
    cofactors: list[str] = Field(default_factory=list)

    # Pathway footprint
    source_kegg_pathways: list[dict] = Field(default_factory=list)   # [{id, name}]
    source_reactome_pathways: list[dict] = Field(default_factory=list)

    # Conflict analysis
    pathway_conflicts: list[PathwayConflict] = Field(default_factory=list)
    competing_host_ec_enzymes: list[dict] = Field(default_factory=list)   # same-EC host enzymes
    metabolite_overlap_pathways: list[dict] = Field(default_factory=list) # host pathways sharing metabolites

    # PPI impact
    host_ppi_partners: list[dict] = Field(default_factory=list)  # STRING predictions in host
    hub_disruption_risk: str = "Unknown"

    # Risk summary (Gemini-generated)
    overall_metabolic_risk: str = ""    # High | Moderate | Low | Unknown
    key_conflicts: list[str] = Field(default_factory=list)
    gemini_analysis: str = ""
    experimental_approaches: list[str] = Field(default_factory=list)
    analysis_caveats: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# KEGG EC-level helpers
# ---------------------------------------------------------------------------

_KEGG_BASE = "https://rest.kegg.jp"


async def _fetch_ec_entry(ec_number: str, client: httpx.AsyncClient) -> dict:
    """Fetch and parse a KEGG EC entry (substrates, products, reaction)."""
    clean = ec_number.strip().lstrip("EC").strip().lstrip(":")
    try:
        r = await client.get(f"{_KEGG_BASE}/get/ec:{clean}", timeout=30)
        if r.status_code != 200 or not r.text.strip():
            return {}
        # Parse flat-file
        result: dict = {}
        current_key = ""
        for line in r.text.split("\n"):
            if line.startswith("///"):
                break
            if line and not line[0].isspace():
                parts = line.split(None, 1)
                current_key = parts[0]
                result[current_key] = parts[1] if len(parts) > 1 else ""
            elif current_key:
                result[current_key] = result[current_key] + " " + line.strip()
        return result
    except Exception:
        return {}


async def _find_human_pathways_for_ec(ec_number: str, client: httpx.AsyncClient) -> list[str]:
    """Find human KEGG pathway IDs that contain enzymes with this EC number."""
    clean = ec_number.strip().lstrip("EC").strip().lstrip(":")
    try:
        r = await client.get(f"{_KEGG_BASE}/link/pathway/ec:{clean}", timeout=30)
        if r.status_code != 200:
            return []
        ids = []
        for line in r.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pid = parts[1].strip()
                if "hsa" in pid:   # human pathways only
                    ids.append(pid)
        return ids
    except Exception:
        return []


async def _get_pathway_name(pathway_id: str, client: httpx.AsyncClient) -> str:
    """Get the human-readable name of a KEGG pathway."""
    clean = pathway_id.replace("path:", "").strip()
    try:
        r = await client.get(f"{_KEGG_BASE}/list/{clean}", timeout=20)
        if r.status_code == 200 and r.text.strip():
            line = r.text.strip().split("\n")[0]
            parts = line.split("\t")
            if len(parts) >= 2:
                return parts[1].strip()
    except Exception:
        pass
    return clean


async def _get_enzymes_in_pathway(pathway_id: str, client: httpx.AsyncClient) -> list[str]:
    """Return EC numbers of all enzymes in a KEGG human pathway."""
    clean = pathway_id.replace("path:", "").strip()
    try:
        r = await client.get(f"{_KEGG_BASE}/link/enzyme/{clean}", timeout=30)
        if r.status_code != 200:
            return []
        ecs = []
        for line in r.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                ec = parts[1].replace("ec:", "").strip()
                if ec:
                    ecs.append(ec)
        return ecs
    except Exception:
        return []


async def _search_compound_in_human_pathways(
    compound_name: str, client: httpx.AsyncClient
) -> list[dict]:
    """Search KEGG COMPOUND by name and find human pathways containing it."""
    results = []
    try:
        r = await client.get(
            f"{_KEGG_BASE}/find/compound/{compound_name.replace(' ', '+')}", timeout=20
        )
        if r.status_code != 200 or not r.text.strip():
            return []
        compound_ids = []
        for line in r.text.strip().split("\n")[:3]:   # top 3 matches
            parts = line.split("\t")
            if parts:
                compound_ids.append(parts[0].replace("cpd:", "").strip())

        for cid in compound_ids:
            pr = await client.get(f"{_KEGG_BASE}/link/pathway/{cid}", timeout=20)
            if pr.status_code == 200:
                for line in pr.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) >= 2 and "hsa" in parts[1]:
                        results.append({"compound": compound_name, "pathway_id": parts[1].strip()})
    except Exception:
        pass
    return results


def _parse_substrates_products(ec_entry: dict) -> tuple[list[str], list[str], list[str]]:
    """Extract substrates, products, and cofactors from a parsed KEGG EC entry."""
    substrates, products, cofactors = [], [], []

    # Known cofactor keywords
    cofactor_kw = {
        "ATP", "ADP", "NAD+", "NADH", "NADP+", "NADPH", "CoA", "FAD", "FADH2",
        "GTP", "FMN", "TPP", "Mg2+", "Mn2+", "Zn2+", "Fe2+", "Fe-S", "Heme",
        "CO2", "O2", "H2O", "Pi", "PPi", "acetyl-CoA", "glutathione",
    }

    all_text = (ec_entry.get("REACTION", "") + " " + ec_entry.get("SUBSTRATE", "") +
                " " + ec_entry.get("PRODUCT", ""))

    # Try SUBSTRATE / PRODUCT sections first
    for token in ec_entry.get("SUBSTRATE", "").split(";"):
        t = token.strip()
        if t:
            if any(kw.lower() in t.lower() for kw in cofactor_kw):
                cofactors.append(t)
            else:
                substrates.append(t)

    for token in ec_entry.get("PRODUCT", "").split(";"):
        t = token.strip()
        if t:
            if any(kw.lower() in t.lower() for kw in cofactor_kw):
                cofactors.append(t)
            else:
                products.append(t)

    # Fall back to REACTION field if sections absent
    if not substrates and not products:
        rxn = ec_entry.get("REACTION", "")
        if "=" in rxn:
            left, right = rxn.split("=", 1)
            for tok in left.split("+"):
                t = tok.strip()
                if t and any(kw.lower() in t.lower() for kw in cofactor_kw):
                    cofactors.append(t)
                elif t:
                    substrates.append(t)
            for tok in right.split("+"):
                t = tok.strip()
                if t and any(kw.lower() in t.lower() for kw in cofactor_kw):
                    cofactors.append(t)
                elif t:
                    products.append(t)

    return (
        list(dict.fromkeys(substrates)),
        list(dict.fromkeys(products)),
        list(dict.fromkeys(cofactors)),
    )


# ---------------------------------------------------------------------------
# STRING interaction lookup (host organism)
# ---------------------------------------------------------------------------

async def _fetch_host_ppi(
    gene_symbol: str,
    target_taxid: int,
    client: httpx.AsyncClient,
    min_score: int = 700,
    limit: int = 15,
) -> list[dict]:
    """Fetch high-confidence STRING interactions for gene_symbol in target organism."""
    try:
        r = await client.get(
            "https://string-db.org/api/json/interaction_partners",
            params={
                "identifiers": gene_symbol,
                "species": target_taxid,
                "limit": limit,
                "required_score": min_score,
                "caller_identity": "proteinscope",
            },
            timeout=30,
        )
        if r.status_code == 200:
            data = r.json()
            return data if isinstance(data, list) else []
    except Exception:
        pass
    return []


# ---------------------------------------------------------------------------
# Reactome pathway lookup
# ---------------------------------------------------------------------------

async def _fetch_source_reactome_pathways(uniprot_id: str, client: httpx.AsyncClient) -> list[dict]:
    """Fetch Reactome pathways for a source UniProt accession."""
    try:
        r = await client.get(
            f"https://reactome.org/ContentService/data/mapping/UniProt/{uniprot_id}/pathways",
            timeout=30,
        )
        if r.status_code == 200:
            data = r.json()
            if isinstance(data, list):
                return [{"id": p.get("stId", ""), "name": p.get("displayName", "")} for p in data[:10]]
    except Exception:
        pass
    return []


# ---------------------------------------------------------------------------
# Conflict classifier
# ---------------------------------------------------------------------------

def _classify_conflicts(
    ec_number: Optional[str],
    substrates: list[str],
    products: list[str],
    cofactors: list[str],
    human_pathway_names: list[str],
    competing_ec_list: list[str],
    source_organelle: Optional[str] = None,
) -> list[PathwayConflict]:
    """Rule-based conflict classifier applied before Gemini synthesis."""
    conflicts: list[PathwayConflict] = []

    # 1. Organelle conflict — missing compartment in target
    if source_organelle == "Chloroplast":
        conflicts.append(PathwayConflict(
            pathway_id="ORGANELLE_CONFLICT",
            pathway_name="Chloroplast-dependent reactions",
            conflict_type="missing_organelle",
            description=(
                "The source protein is localized to the chloroplast, an organelle absent "
                "in animal cells. Without a chloroplast, the protein cannot assemble into "
                "its native complex, lacks its electron transport partner proteins, and the "
                "substrate ribulose-1,5-bisphosphate (RuBP) is not synthesized in human cells."
            ),
            severity="High",
            evidence_sources=["UniProt subcellular location"],
        ))

    # 2. EC number duplication — direct catalytic competition in host
    if ec_number and competing_ec_list:
        conflicts.append(PathwayConflict(
            pathway_id="EC_COMPETITION",
            pathway_name="Competing enzymatic reactions",
            conflict_type="substrate_competition",
            host_genes_affected=competing_ec_list[:5],
            description=(
                f"Human enzymes with the same EC number ({ec_number}) already catalyze "
                f"this reaction class. Ectopic expression of the foreign enzyme would "
                f"compete for substrate flux and could dysregulate existing metabolic homeostasis. "
                f"Affected host enzymes: {', '.join(competing_ec_list[:5])}."
            ),
            severity="High",
            evidence_sources=["KEGG EC linkage"],
        ))

    # 3. Cofactor depletion — critical cofactors that are shared across many host pathways
    critical_cofactors = {"NAD+", "NADH", "NADP+", "NADPH", "ATP", "ADP", "CoA", "acetyl-CoA"}
    depleted = [c for c in cofactors if any(kc in c for kc in critical_cofactors)]
    if depleted:
        conflicts.append(PathwayConflict(
            pathway_id="COFACTOR_DEPLETION",
            pathway_name="Cofactor-sharing host pathways",
            conflict_type="cofactor_depletion",
            description=(
                f"The inserted enzyme requires {', '.join(depleted[:3])} as cofactor(s). "
                f"These are rate-limiting molecules shared across hundreds of host metabolic "
                f"reactions. Constitutive high-level expression of the insert could create "
                f"cofactor competition and deplete availability for essential host reactions."
            ),
            severity="Moderate",
            evidence_sources=["KEGG EC entry"],
        ))

    # 4. CO2-fixing enzymes — specific conflict with human carbonic anhydrase system
    if ec_number and ec_number.startswith("4.1.1."):
        conflicts.append(PathwayConflict(
            pathway_id="CO2_BALANCE",
            pathway_name="Carbonic anhydrase / CO₂-HCO₃⁻ equilibrium",
            conflict_type="substrate_competition",
            host_genes_affected=["CA1", "CA2", "CA4", "CA9"],
            description=(
                "Carboxylase activity (EC 4.1.1.x) consumes CO₂. In human cells, CO₂ "
                "is a tightly regulated molecule controlling blood pH via carbonic anhydrase "
                "(CA1/CA2). Ectopic CO₂ fixation would disrupt the CO₂-HCO₃⁻-H⁺ equilibrium, "
                "potentially causing intracellular acidification."
            ),
            severity="Moderate",
            evidence_sources=["KEGG reaction chemistry"],
        ))

    # 5. Products appearing in human pathways → feedback / accumulation risk
    glycolysis_metabolites = {
        "3-phosphoglycerate", "3-PGA", "pyruvate", "fructose-6-phosphate",
        "glucose-6-phosphate", "phosphoenolpyruvate", "2-phosphoglycerate",
    }
    for prod in products:
        for gly_met in glycolysis_metabolites:
            if gly_met.lower() in prod.lower():
                conflicts.append(PathwayConflict(
                    pathway_id="GLYCOLYSIS_INTERFERENCE",
                    pathway_name="Glycolysis / Gluconeogenesis (hsa00010)",
                    conflict_type="product_inhibition",
                    description=(
                        f"The inserted enzyme produces {prod}, a metabolite that already "
                        f"participates in human glycolysis/gluconeogenesis. Excess production "
                        f"could allosterically inhibit upstream glycolytic enzymes (product "
                        f"feedback inhibition) or drive futile cycles."
                    ),
                    severity="Moderate",
                    evidence_sources=["KEGG compound mapping"],
                ))
                break

    # 6. Oxygenase side-reaction toxicity (RuBisCO-type)
    if ec_number in ("4.1.1.39", "1.1.1.39"):
        conflicts.append(PathwayConflict(
            pathway_id="PHOTORESPIRATION_TOXIN",
            pathway_name="Photorespiratory pathway (absent in human cells)",
            conflict_type="product_toxicity",
            host_genes_affected=["SHMT1", "SHMT2", "AGXT"],
            description=(
                "RuBisCO has both carboxylase and oxygenase activity. In the presence of O₂ "
                "(always present in human cells), it produces 2-phosphoglycolate — a toxic "
                "metabolite that inhibits triose-phosphate isomerase (TPI) and disrupts "
                "glycolytic flux. In plants this is detoxified by photorespiration; "
                "no equivalent pathway exists in human epithelial cells."
            ),
            severity="High",
            evidence_sources=["Biochemical literature", "KEGG reaction R00024"],
        ))

    return conflicts


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_impact_with_gemini(report: InsertionImpactReport) -> InsertionImpactReport:
    """Use Gemini 2.5 Pro to synthesise a scientific conflict analysis from all computed facts."""
    try:
        from core.gemini_interpreter import _call

        # Build fact block
        lines: list[str] = [
            f"Gene being inserted: {report.source_gene} from {report.source_organism}",
            f"Host organism: {report.target_organism}",
        ]
        if report.source_ec_number:
            lines.append(f"Enzyme Commission number: {report.source_ec_number}")
        if report.reaction_equation:
            lines.append(f"Catalyzed reaction: {report.reaction_equation}")
        if report.substrates:
            lines.append(f"Substrates: {', '.join(report.substrates[:6])}")
        if report.products:
            lines.append(f"Products: {', '.join(report.products[:6])}")
        if report.cofactors:
            lines.append(f"Required cofactors: {', '.join(report.cofactors[:6])}")
        if report.source_kegg_pathways:
            names = [p.get("name", p.get("id", "?")) for p in report.source_kegg_pathways[:5]]
            lines.append(f"Source organism pathways: {'; '.join(names)}")
        if report.competing_host_ec_enzymes:
            ecs = [e.get("ec", "?") for e in report.competing_host_ec_enzymes[:5]]
            lines.append(f"Host enzymes with same EC number: {', '.join(ecs)}")
        if report.metabolite_overlap_pathways:
            mop = [p.get("pathway_id", "?") for p in report.metabolite_overlap_pathways[:6]]
            lines.append(f"Host pathways sharing metabolites with insert: {', '.join(mop)}")
        if report.host_ppi_partners:
            partners = [p.get("preferredName_B", p.get("stringId_B", "?"))
                       for p in report.host_ppi_partners[:8]]
            lines.append(f"Predicted host PPI partners (STRING): {', '.join(partners)}")
        if report.pathway_conflicts:
            lines.append("\nPre-computed conflicts:")
            for c in report.pathway_conflicts:
                lines.append(f"  [{c.severity}] {c.conflict_type}: {c.description[:150]}…")
        if report.source_reactome_pathways:
            rnames = [p.get("name", "?") for p in report.source_reactome_pathways[:4]]
            lines.append(f"Reactome pathways: {'; '.join(rnames)}")

        facts_block = "\n".join(lines)

        prompt = (
            "You are a senior systems biologist and metabolic engineer. "
            "Analyze the following data about inserting a foreign gene into a host organism "
            "and produce a rigorous, evidence-based metabolic impact assessment.\n\n"
            f"{facts_block}\n\n"
            "Using standard bioinformatics analysis methods, assess:\n"
            "  1. METABOLIC PATHWAY CONFLICTS — which specific host pathways are disrupted and how?\n"
            "  2. ENZYME INHIBITION — which host enzymes are inhibited by substrates/products/cofactors of the insert?\n"
            "  3. NETWORK PERTURBATION — how does PPI disruption propagate through the host metabolic network?\n"
            "  4. THERMODYNAMIC FEASIBILITY — would the inserted reaction actually proceed in host cell conditions?\n"
            "  5. OVERALL RISK LEVEL — High / Moderate / Low with justification\n\n"
            "Return ONLY a raw JSON object (no markdown fences):\n"
            "{\n"
            "  \"overall_metabolic_risk\": \"High|Moderate|Low|Unknown\",\n"
            "  \"key_conflicts\": [\"<concise conflict 1>\", \"<concise conflict 2>\", ...],\n"
            "  \"gemini_analysis\": \"<5-8 sentences scientifically rigorous analysis referencing the evidence>\",\n"
            "  \"experimental_approaches\": [\"<specific experiment 1>\", ...],\n"
            "  \"analysis_caveats\": [\"<caveat 1>\", ...]\n"
            "}"
        )

        raw = await _call(prompt)
        if not raw:
            report.overall_metabolic_risk = "Unknown"
            report.gemini_analysis = "AI synthesis unavailable."
            return report

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        valid_risk = {"High", "Moderate", "Low", "Unknown"}
        risk = str(data.get("overall_metabolic_risk", "Unknown")).strip()
        report.overall_metabolic_risk = risk if risk in valid_risk else "Unknown"

        kc = data.get("key_conflicts", [])
        report.key_conflicts = [str(x).strip() for x in kc if x] if isinstance(kc, list) else []

        report.gemini_analysis = str(data.get("gemini_analysis", "")).strip()

        ea = data.get("experimental_approaches", [])
        report.experimental_approaches = [str(x).strip() for x in ea if x] if isinstance(ea, list) else []

        ac = data.get("analysis_caveats", [])
        report.analysis_caveats = [str(x).strip() for x in ac if x] if isinstance(ac, list) else []

    except Exception:
        if not report.overall_metabolic_risk:
            report.overall_metabolic_risk = "Unknown"
            report.gemini_analysis = "AI synthesis failed."

    return report


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

_HUMAN_TAXID = 9606


async def run_insertion_impact_analysis(
    source_gene: str,
    source_organism: str,
    target_organism: str,
    source_uniprot_id: Optional[str] = None,
    source_ec_number: Optional[str] = None,
    source_organelle: Optional[str] = None,
    step_cb=None,
) -> InsertionImpactReport:
    """Analyze the metabolic impact of inserting source_gene into target_organism.

    Args:
        source_gene:        Gene symbol (e.g. "RBCL").
        source_organism:    Donor organism (e.g. "plant").
        target_organism:    Host organism (e.g. "Homo sapiens").
        source_uniprot_id:  Optional UniProt accession of source protein.
        source_ec_number:   Optional pre-resolved EC number.
        source_organelle:   Optional organelle context (e.g. "Chloroplast").
        step_cb:            Optional async progress callback.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    report = InsertionImpactReport(
        source_gene=source_gene.upper(),
        source_organism=source_organism,
        target_organism=target_organism,
        source_ec_number=source_ec_number,
        source_organelle=source_organelle,
    )

    async with httpx.AsyncClient(timeout=30) as client:

        # ── 1. Resolve EC number if not provided ────────────────────────────
        await _step("[1/7] Resolving enzyme EC number…")
        if not report.source_ec_number and source_uniprot_id:
            try:
                from fetchers.brenda import get_ec_number
                from fetchers import uniprot as uni
                entry = await uni.fetch_by_accession(source_uniprot_id)
                if isinstance(entry, dict):
                    report.source_ec_number = get_ec_number(entry)
            except Exception:
                pass

        # Fallback: try MetaCyc
        if not report.source_ec_number:
            try:
                from fetchers.metacyc import fetch_enzyme_reactions, _parse_ec_numbers
                rxn_xml = await fetch_enzyme_reactions(source_gene)
                ec_list = _parse_ec_numbers(rxn_xml)
                if ec_list:
                    report.source_ec_number = ec_list[0]
            except Exception:
                pass

        # ── 2. KEGG EC entry — substrates, products, reaction ───────────────
        await _step("[2/7] Fetching KEGG enzyme entry (substrates, products, cofactors)…")
        ec_entry: dict = {}
        if report.source_ec_number:
            ec_entry = await _fetch_ec_entry(report.source_ec_number, client)
            if ec_entry:
                report.reaction_equation = ec_entry.get("REACTION", "").strip() or None
                subs, prods, cofs = _parse_substrates_products(ec_entry)
                report.substrates = subs
                report.products = prods
                report.cofactors = cofs

        # ── 3. Human KEGG pathways with same EC → competing enzymes ─────────
        await _step("[3/7] Mapping competing host enzymes via KEGG pathway/EC linkage…")
        if report.source_ec_number:
            try:
                human_pathway_ids = await _find_human_pathways_for_ec(
                    report.source_ec_number, client
                )
                # Get names and competing enzyme EC lists
                comp_ec_set: set[str] = set()
                for pid in human_pathway_ids[:5]:
                    name = await _get_pathway_name(pid, client)
                    report.source_kegg_pathways.append({"id": pid, "name": name})
                    ecs = await _get_enzymes_in_pathway(pid, client)
                    for ec in ecs:
                        if ec != report.source_ec_number:
                            comp_ec_set.add(ec)
                report.competing_host_ec_enzymes = [{"ec": ec} for ec in list(comp_ec_set)[:10]]
            except Exception:
                pass

        # ── 4. Metabolite overlap — host pathways sharing substrates/products ─
        await _step("[4/7] Searching host pathways for metabolite overlap…")
        overlap_set: dict[str, dict] = {}
        all_metabolites = (report.substrates + report.products)[:6]
        overlap_tasks = [
            _search_compound_in_human_pathways(met, client)
            for met in all_metabolites
        ]
        try:
            overlap_results = await asyncio.gather(*overlap_tasks, return_exceptions=True)
            for result in overlap_results:
                if isinstance(result, list):
                    for item in result:
                        pid = item.get("pathway_id", "")
                        if pid and pid not in overlap_set:
                            name = await _get_pathway_name(pid, client)
                            overlap_set[pid] = {"pathway_id": pid, "name": name,
                                               "compound": item.get("compound", "")}
            report.metabolite_overlap_pathways = list(overlap_set.values())[:8]
        except Exception:
            pass

        # ── 5. STRING — host PPI predictions ────────────────────────────────
        await _step("[5/7] Predicting host protein-protein interactions (STRING)…")
        try:
            ppi = await _fetch_host_ppi(source_gene, _HUMAN_TAXID, client, min_score=700)
            report.host_ppi_partners = ppi[:12]
            if ppi:
                # Simple hub detection: if >5 high-confidence partners, flag as hub risk
                report.hub_disruption_risk = "High" if len(ppi) >= 5 else "Moderate"
        except Exception:
            pass

        # ── 6. Reactome — source protein pathway context ─────────────────────
        await _step("[6/7] Fetching Reactome pathway context…")
        if source_uniprot_id:
            try:
                report.source_reactome_pathways = await _fetch_source_reactome_pathways(
                    source_uniprot_id, client
                )
            except Exception:
                pass

        # ── 7. Rule-based conflict classification ────────────────────────────
        await _step("[7/7] Classifying metabolic conflicts…")
        human_pathway_names = [p.get("name", "") for p in report.source_kegg_pathways]
        competing_ecs = [e.get("ec", "") for e in report.competing_host_ec_enzymes]
        report.pathway_conflicts = _classify_conflicts(
            ec_number=report.source_ec_number,
            substrates=report.substrates,
            products=report.products,
            cofactors=report.cofactors,
            human_pathway_names=human_pathway_names,
            competing_ec_list=competing_ecs,
            source_organelle=source_organelle,
        )

    # ── Gemini synthesis (outside the httpx context) ─────────────────────
    await _step("[7/7] Gemini 2.5 Pro synthesizing conflict analysis…")
    report = await _synthesize_impact_with_gemini(report)

    return report
