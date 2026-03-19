"""Reverse Genetics Analyzer.

Given a gene and mutation type, predicts the downstream metabolic and pathological
consequences by cross-referencing all connected databases.

Analysis pipeline:
1.  Resolve gene → UniProt entry (protein function, domains, active sites, EC number)
2.  Classify mutation impact on protein structure/function (domain-aware)
3.  ClinVar  → known pathogenic variants, clinical significance, associated conditions
4.  OMIM     → mendelian disease associations, inheritance patterns
5.  KEGG     → disrupted metabolic pathways, compensatory routes
6.  Reactome → disrupted biological processes and cascades
7.  STRING   → lost protein-protein interactions, network centrality impact
8.  PharmGKB → altered drug response (PGx consequences of mutation)
9.  NCBI Gene → paralog/compensatory gene identification
10. Gemini 2.5 Pro → integrative bioinformatics synthesis

Mutation types handled:
  loss_of_function (LoF)  — knockout, frameshift, nonsense, splice-site, deletion
  gain_of_function (GoF)  — activating missense, amplification, overexpression
  missense                — single amino acid substitution (domain-aware impact)
  dominant_negative       — mutant sequesters WT complexes
  haploinsufficiency      — one functional copy is insufficient
  general                 — unspecified; Gemini infers likely mechanism
"""

from __future__ import annotations

import asyncio
import json as _json
from typing import Optional

import httpx
from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class PathwayConsequence(BaseModel):
    pathway_id: str
    pathway_name: str
    consequence_type: str   # disrupted | hyperactivated | compensated | bypassed
    affected_reactions: list[str] = Field(default_factory=list)
    downstream_effects: list[str] = Field(default_factory=list)
    severity: str           # High | Moderate | Low
    source: str             # KEGG | Reactome | STRING


class DomainImpact(BaseModel):
    domain_name: str
    domain_type: str        # active_site | binding_site | signal | transmembrane | pfam_domain
    mutation_in_domain: bool
    functional_consequence: str


class ClinicalCorrelation(BaseModel):
    disease: str
    significance: str       # Pathogenic | Likely pathogenic | VUS | Benign
    inheritance: Optional[str] = None   # AD | AR | XL | mitochondrial
    omim_id: Optional[str] = None
    source: str             # ClinVar | OMIM | PharmGKB


class ReverseGeneticsReport(BaseModel):
    gene: str
    organism: str
    mutation_type: str
    specific_variant: Optional[str] = None   # e.g. "p.R175H", "c.185delAG"
    target_domain: Optional[str] = None      # if variant maps to a known domain

    # Protein context
    protein_name: Optional[str] = None
    protein_function: Optional[str] = None
    ec_number: Optional[str] = None
    uniprot_id: Optional[str] = None
    is_essential: bool = False              # inferred from essentiality databases
    is_tumor_suppressor: bool = False
    is_oncogene: bool = False

    # Domain impact analysis
    domain_impacts: list[DomainImpact] = Field(default_factory=list)

    # Pathway consequences (KEGG + Reactome)
    disrupted_kegg_pathways: list[dict] = Field(default_factory=list)    # [{id, name, role}]
    disrupted_reactome_pathways: list[dict] = Field(default_factory=list)
    pathway_consequences: list[PathwayConsequence] = Field(default_factory=list)

    # PPI network impact (STRING)
    lost_interactions: list[dict] = Field(default_factory=list)  # high-conf partners lost
    hub_score: Optional[float] = None          # fraction of network affected
    hub_classification: str = "Unknown"        # Essential hub | Peripheral | Unknown

    # Clinical associations
    known_clinvar_variants: list[dict] = Field(default_factory=list)
    clinical_correlations: list[ClinicalCorrelation] = Field(default_factory=list)
    omim_diseases: list[str] = Field(default_factory=list)

    # PGx implications
    pgx_drug_interactions: list[dict] = Field(default_factory=list)
    pgx_variants: list[dict] = Field(default_factory=list)

    # Compensatory mechanisms
    paralogs: list[str] = Field(default_factory=list)
    compensatory_pathways: list[str] = Field(default_factory=list)

    # Gemini synthesis
    overall_impact: str = ""            # High | Moderate | Low | Unknown
    key_consequences: list[str] = Field(default_factory=list)
    gemini_analysis: str = ""
    experimental_approaches: list[str] = Field(default_factory=list)
    analysis_caveats: list[str] = Field(default_factory=list)

    # Foundational references
    references: list[dict] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Curated foundational references for Reverse Genetics
# ---------------------------------------------------------------------------

_RG_REFERENCES: list[dict] = [
    {
        "pubmed_id": "22745249",
        "doi": "10.1126/science.1225829",
        "title": "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity",
        "authors": ["Jinek M", "Chylinski K", "Fonfara I", "Hauer M", "Doudna JA", "Charpentier E"],
        "journal": "Science",
        "year": 2012,
        "apa_citation": "Jinek, M., Chylinski, K., Fonfara, I., Hauer, M., Doudna, J. A., & Charpentier, E. (2012). A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. Science, 337(6096), 816–821. https://doi.org/10.1126/science.1225829",
    },
    {
        "pubmed_id": "23287718",
        "doi": "10.1126/science.1231143",
        "title": "Multiplex genome engineering using CRISPR/Cas systems",
        "authors": ["Cong L", "Ran FA", "Cox D", "Lin S", "Barretto R", "Habib N", "Zhang F"],
        "journal": "Science",
        "year": 2013,
        "apa_citation": "Cong, L., Ran, F. A., Cox, D., Lin, S., Barretto, R., Habib, N., & Zhang, F. (2013). Multiplex genome engineering using CRISPR/Cas systems. Science, 339(6121), 819–823. https://doi.org/10.1126/science.1231143",
    },
    {
        "pubmed_id": "24157548",
        "doi": "10.1038/nprot.2013.143",
        "title": "Genome engineering using the CRISPR-Cas9 system",
        "authors": ["Ran FA", "Hsu PD", "Wright J", "Agarwala V", "Scott DA", "Zhang F"],
        "journal": "Nature Protocols",
        "year": 2013,
        "apa_citation": "Ran, F. A., Hsu, P. D., Wright, J., Agarwala, V., Scott, D. A., & Zhang, F. (2013). Genome engineering using the CRISPR-Cas9 system. Nature Protocols, 8(11), 2281–2308. https://doi.org/10.1038/nprot.2013.143",
    },
    {
        "pubmed_id": "31634902",
        "doi": "10.1038/s41586-019-1711-4",
        "title": "Search-and-replace genome editing without double-strand breaks or donor DNA",
        "authors": ["Anzalone AV", "Randolph PB", "Davis JR", "Sousa AA", "Koblan LW", "Levy JM", "Liu DR"],
        "journal": "Nature",
        "year": 2019,
        "apa_citation": "Anzalone, A. V., Randolph, P. B., Davis, J. R., Sousa, A. A., Koblan, L. W., Levy, J. M., & Liu, D. R. (2019). Search-and-replace genome editing without double-strand breaks or donor DNA. Nature, 576(7785), 149–157. https://doi.org/10.1038/s41586-019-1711-4",
    },
    {
        "pubmed_id": "36408920",
        "doi": "10.1093/nar/gkac1052",
        "title": "UniProt: the Universal Protein Knowledgebase in 2023",
        "authors": ["UniProt Consortium"],
        "journal": "Nucleic Acids Research",
        "year": 2023,
        "apa_citation": "UniProt Consortium. (2023). UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052",
    },
]


# ---------------------------------------------------------------------------
# Domain / feature impact analysis
# ---------------------------------------------------------------------------

_FUNCTIONAL_FEATURE_TYPES = {
    "active site", "binding site", "metal binding", "site",
    "disulfide bond", "modified residue",
}

_STRUCTURAL_FEATURE_TYPES = {
    "transmembrane", "signal peptide", "transit peptide",
    "propeptide", "region",
}

_DOMAIN_FEATURE_TYPES = {"domain", "repeat", "zinc finger"}


def _analyze_domain_impacts(
    uniprot_entry: dict,
    specific_variant: Optional[str],
    mutation_type: str,
) -> tuple[list[DomainImpact], Optional[str]]:
    """Infer which domains are impacted by the mutation.

    Returns (list of DomainImpact, target_domain_name or None).
    If specific_variant is None, assumes the mutation type applies globally.
    """
    impacts: list[DomainImpact] = []
    target_domain: Optional[str] = None

    # Parse variant position (e.g. p.R175H → 175, p.del → None)
    variant_pos: Optional[int] = None
    if specific_variant:
        import re
        m = re.search(r"[A-Za-z*]?(\d+)[A-Za-z*]?$", specific_variant.replace("p.", "").replace("c.", ""))
        if m:
            try:
                variant_pos = int(m.group(1))
            except ValueError:
                pass

    features = uniprot_entry.get("features", [])
    for feat in features:
        ftype = feat.get("type", "").lower()
        fname = feat.get("description", ftype) or ftype
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")

        try:
            fstart = int(start) if start is not None else None
            fend = int(end) if end is not None else None
        except (TypeError, ValueError):
            fstart = fend = None

        # Determine if variant falls in this feature
        in_feature = False
        if variant_pos and fstart is not None and fend is not None:
            in_feature = fstart <= variant_pos <= fend
        elif not variant_pos:
            # LoF/GoF — globally affects all functional features
            in_feature = ftype in _FUNCTIONAL_FEATURE_TYPES or ftype in _DOMAIN_FEATURE_TYPES

        if ftype in _FUNCTIONAL_FEATURE_TYPES and (in_feature or not variant_pos):
            if ftype == "active site":
                consequence = (
                    "Direct abolition of catalytic activity" if in_feature
                    else "Possible active site disruption via conformational change"
                )
                dtype = "active_site"
            elif "binding" in ftype:
                consequence = "Loss of substrate/cofactor binding" if in_feature else "Binding site disruption risk"
                dtype = "binding_site"
            else:
                consequence = "Functional modification site disrupted" if in_feature else "Post-translational modification affected"
                dtype = "active_site"

            impacts.append(DomainImpact(
                domain_name=fname,
                domain_type=dtype,
                mutation_in_domain=in_feature,
                functional_consequence=consequence,
            ))
            if in_feature and not target_domain:
                target_domain = fname

        elif ftype in _DOMAIN_FEATURE_TYPES and in_feature:
            impacts.append(DomainImpact(
                domain_name=fname,
                domain_type="pfam_domain",
                mutation_in_domain=True,
                functional_consequence=f"Structural disruption of {fname} domain — folding and interaction interfaces at risk",
            ))
            if not target_domain:
                target_domain = fname

    return impacts, target_domain


def _infer_gene_role(uniprot_entry: dict) -> tuple[bool, bool, bool]:
    """Infer essential / tumor_suppressor / oncogene from UniProt keywords."""
    keywords = [kw.get("name", "").lower() for kw in uniprot_entry.get("keywords", [])]
    comments = " ".join(
        " ".join(t.get("value", "") for t in c.get("texts", []))
        for c in uniprot_entry.get("comments", [])
    ).lower()

    is_essential = any(k in keywords for k in ["essential", "viability"])
    is_tumor_suppressor = any(k in keywords for k in ["tumor suppressor"]) or "tumor suppressor" in comments
    is_oncogene = any(k in keywords for k in ["proto-oncogene"]) or "proto-oncogene" in comments or "oncogene" in comments

    return is_essential, is_tumor_suppressor, is_oncogene


def _extract_protein_function(uniprot_entry: dict) -> str:
    """Get the Function comment from a UniProt entry."""
    for comment in uniprot_entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value", "")
    return ""


# ---------------------------------------------------------------------------
# KEGG pathway consequence inference
# ---------------------------------------------------------------------------

async def _fetch_kegg_pathways_for_gene(ncbi_gene_id: str, client: httpx.AsyncClient) -> list[dict]:
    """Return [{id, name}] for all KEGG human pathways containing this gene."""
    results = []
    try:
        r = await client.get(
            f"https://rest.kegg.jp/link/pathway/hsa:{ncbi_gene_id}", timeout=30
        )
        if r.status_code != 200:
            return []
        pathway_ids = [
            line.split("\t")[1].strip()
            for line in r.text.strip().split("\n")
            if "\t" in line
        ]
        for pid in pathway_ids[:8]:
            r2 = await client.get(f"https://rest.kegg.jp/list/{pid.replace('path:','')}", timeout=20)
            name = pid
            if r2.status_code == 200 and r2.text.strip():
                parts = r2.text.strip().split("\n")[0].split("\t")
                if len(parts) >= 2:
                    name = parts[1].strip()
            results.append({"id": pid, "name": name})
    except Exception:
        pass
    return results


async def _build_pathway_consequences(
    pathways: list[dict],
    mutation_type: str,
    gene: str,
) -> list[PathwayConsequence]:
    """Convert pathway list into structured pathway consequence objects."""
    lof_types = {"loss_of_function", "dominant_negative", "haploinsufficiency", "frameshift", "nonsense", "deletion"}
    gof_types = {"gain_of_function", "amplification"}

    consequences = []
    for pw in pathways:
        pid = pw.get("id", "")
        pname = pw.get("name", "?")

        if mutation_type.lower() in lof_types:
            ctype = "disrupted"
            severity = "High"
            downstream = [
                f"Loss of {gene}'s contribution to {pname}",
                "Pathway flux reduction or block at reaction step",
                "Potential accumulation of upstream metabolites",
                "Deficiency of downstream products/signals",
            ]
        elif mutation_type.lower() in gof_types:
            ctype = "hyperactivated"
            severity = "Moderate"
            downstream = [
                f"Excess {gene} activity drives {pname} beyond normal flux",
                "Risk of downstream product accumulation and toxicity",
                "Negative feedback circuits may be overwhelmed",
            ]
        else:
            ctype = "potentially disrupted"
            severity = "Moderate"
            downstream = [
                f"Partial loss of {gene} function in {pname}",
                "Severity depends on mutation position relative to active site",
            ]

        consequences.append(PathwayConsequence(
            pathway_id=pid,
            pathway_name=pname,
            consequence_type=ctype,
            affected_reactions=[],
            downstream_effects=downstream,
            severity=severity,
            source="KEGG",
        ))

    return consequences


# ---------------------------------------------------------------------------
# STRING network impact
# ---------------------------------------------------------------------------

async def _fetch_ppi_partners(gene: str, taxid: int, client: httpx.AsyncClient) -> list[dict]:
    """Fetch high-confidence STRING interaction partners."""
    try:
        r = await client.get(
            "https://string-db.org/api/json/interaction_partners",
            params={
                "identifiers": gene,
                "species": taxid,
                "limit": 20,
                "required_score": 700,
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


def _classify_hub(interaction_count: int) -> tuple[str, float]:
    """Classify hub importance from interaction count."""
    if interaction_count >= 15:
        return "Essential hub — loss will broadly disrupt PPI network", min(1.0, interaction_count / 20)
    elif interaction_count >= 5:
        return "Moderate hub — several pathways will be affected", interaction_count / 20
    elif interaction_count > 0:
        return "Peripheral node — limited network impact", interaction_count / 20
    return "Unknown — no high-confidence interactions found", 0.0


# ---------------------------------------------------------------------------
# Paralog/compensatory gene search (NCBI gene search)
# ---------------------------------------------------------------------------

async def _find_paralogs(gene: str, client: httpx.AsyncClient) -> list[str]:
    """Try to find paralogs via NCBI gene homolog endpoint."""
    try:
        r = await client.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": "gene", "term": f"{gene}[gene] AND 9606[taxid]", "retmode": "json"},
            timeout=20,
        )
        if r.status_code != 200:
            return []
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        gene_id = ids[0]

        r2 = await client.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
            params={
                "dbfrom": "gene",
                "db": "gene",
                "id": gene_id,
                "linkname": "gene_gene_homologs",
                "retmode": "json",
            },
            timeout=20,
        )
        if r2.status_code != 200:
            return []
        links = r2.json().get("linksets", [{}])[0].get("linksetdbs", [{}])[0].get("links", [])
        # Fetch gene symbols for paralog IDs (limit 5)
        paralogs = []
        for lid in links[:5]:
            r3 = await client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={"db": "gene", "id": lid, "retmode": "json"},
                timeout=15,
            )
            if r3.status_code == 200:
                try:
                    sym = r3.json().get("result", {}).get(str(lid), {}).get("name", "")
                    if sym and sym != gene:
                        paralogs.append(sym)
                except Exception:
                    pass
        return paralogs
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Gemini synthesis — bioinformatics-level analysis
# ---------------------------------------------------------------------------

async def _synthesize_with_gemini(report: ReverseGeneticsReport) -> ReverseGeneticsReport:
    try:
        from core.gemini_interpreter import _call

        lines: list[str] = [
            f"Gene: {report.gene} ({report.organism})",
            f"Mutation type: {report.mutation_type}",
        ]
        if report.specific_variant:
            lines.append(f"Specific variant: {report.specific_variant}")
        if report.protein_name:
            lines.append(f"Protein: {report.protein_name}")
        if report.protein_function:
            lines.append(f"Normal function: {report.protein_function[:300]}")
        if report.ec_number:
            lines.append(f"EC number: {report.ec_number}")

        role_flags = []
        if report.is_essential:
            role_flags.append("ESSENTIAL GENE")
        if report.is_tumor_suppressor:
            role_flags.append("TUMOR SUPPRESSOR")
        if report.is_oncogene:
            role_flags.append("PROTO-ONCOGENE")
        if role_flags:
            lines.append(f"Gene role: {', '.join(role_flags)}")

        if report.domain_impacts:
            affected = [d for d in report.domain_impacts if d.mutation_in_domain]
            if affected:
                lines.append("DIRECT DOMAIN IMPACTS:")
                for d in affected[:4]:
                    lines.append(f"  - [{d.domain_type}] {d.domain_name}: {d.functional_consequence}")

        if report.disrupted_kegg_pathways:
            names = [p.get("name", "?") for p in report.disrupted_kegg_pathways[:6]]
            lines.append(f"Disrupted KEGG pathways ({len(report.disrupted_kegg_pathways)} total): {'; '.join(names)}")

        if report.disrupted_reactome_pathways:
            names = [p.get("name", "?") for p in report.disrupted_reactome_pathways[:4]]
            lines.append(f"Disrupted Reactome pathways: {'; '.join(names)}")

        if report.lost_interactions:
            partners = [p.get("preferredName_B", "?") for p in report.lost_interactions[:8]]
            lines.append(f"Lost PPI partners (STRING ≥ 0.7): {', '.join(partners)}")
            lines.append(f"Hub classification: {report.hub_classification}")

        if report.clinical_correlations:
            diseases = list({c.disease for c in report.clinical_correlations[:6]})
            lines.append(f"Known clinical associations: {'; '.join(diseases)}")

        if report.omim_diseases:
            lines.append(f"OMIM disease associations: {'; '.join(report.omim_diseases[:4])}")

        if report.pgx_drug_interactions:
            drugs = [d.get("drug_name", d.get("name", "?")) for d in report.pgx_drug_interactions[:5]]
            lines.append(f"Affected drug responses (PharmGKB): {', '.join(drugs)}")

        if report.paralogs:
            lines.append(f"Potential compensatory paralogs: {', '.join(report.paralogs)}")

        facts_block = "\n".join(lines)

        mut_guidance = {
            "loss_of_function": (
                "Address: which pathway steps are ablated, whether paralogs can compensate, "
                "what metabolites accumulate upstream, what downstream signals are lost, "
                "and which diseases arise from LoF in this gene."
            ),
            "gain_of_function": (
                "Address: which pathways are constitutively activated, what feedback loops are bypassed, "
                "oncogenic consequences, toxic metabolite accumulation, and therapeutic implications."
            ),
            "missense": (
                "Address: whether the mutated residue is in a functional domain, how this affects "
                "protein stability and activity, partial vs. complete LoF, dominant-negative potential, "
                "and genotype-phenotype correlations from ClinVar."
            ),
            "dominant_negative": (
                "Address: how the mutant protein sequesters WT complexes, which interaction partners "
                "are trapped, the effective LoF of wild-type copies, and haploinsufficiency phenotypes."
            ),
        }
        specific_guidance = mut_guidance.get(
            report.mutation_type.lower(),
            "Infer the most likely mechanism and address pathway, network, and clinical consequences."
        )

        prompt = (
            "You are a senior molecular geneticist and systems biologist conducting "
            "a reverse genetics analysis. Using ALL the following database-derived evidence, "
            "produce a rigorous, evidence-grounded bioinformatics analysis.\n\n"
            f"EVIDENCE:\n{facts_block}\n\n"
            f"ANALYSIS REQUIREMENTS: {specific_guidance}\n\n"
            "Structure your analysis using standard bioinformatics reasoning:\n"
            "  1. Molecular mechanism of the mutation (what breaks at the protein level)\n"
            "  2. Primary pathway consequences (direct metabolic/signaling disruption)\n"
            "  3. Secondary/cascade effects (indirect consequences through the network)\n"
            "  4. Clinical/pathological phenotype prediction\n"
            "  5. Compensatory mechanisms (paralogs, alternative pathways)\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            "{\n"
            "  \"overall_impact\": \"High|Moderate|Low|Unknown\",\n"
            "  \"key_consequences\": [\"<consequence 1>\", \"<consequence 2>\", ...],\n"
            "  \"gemini_analysis\": \"<6-10 sentence scientifically rigorous analysis>\",\n"
            "  \"experimental_approaches\": [\"<specific experiment 1>\", ...],\n"
            "  \"analysis_caveats\": [\"<caveat 1>\", ...]\n"
            "}"
        )

        raw = await _call(prompt)
        if not raw:
            report.overall_impact = "Unknown"
            report.gemini_analysis = "AI synthesis unavailable."
            return report

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        valid = {"High", "Moderate", "Low", "Unknown"}
        impact = str(data.get("overall_impact", "Unknown")).strip()
        report.overall_impact = impact if impact in valid else "Unknown"

        kc = data.get("key_consequences", [])
        report.key_consequences = [str(x).strip() for x in kc if x] if isinstance(kc, list) else []

        report.gemini_analysis = str(data.get("gemini_analysis", "")).strip()

        ea = data.get("experimental_approaches", [])
        report.experimental_approaches = [str(x).strip() for x in ea if x] if isinstance(ea, list) else []

        ac = data.get("analysis_caveats", [])
        report.analysis_caveats = [str(x).strip() for x in ac if x] if isinstance(ac, list) else []

    except Exception:
        if not report.overall_impact:
            report.overall_impact = "Unknown"
            report.gemini_analysis = "AI synthesis failed due to an error."

    return report


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

_HUMAN_TAXID = 9606


async def run_reverse_genetics_analysis(
    gene: str,
    mutation_type: str = "loss_of_function",
    specific_variant: Optional[str] = None,
    organism: str = "Homo sapiens",
    step_cb=None,
) -> ReverseGeneticsReport:
    """Run full reverse genetics analysis for a gene mutation.

    Args:
        gene:             Gene symbol (e.g. "BRCA1", "TP53", "POLD1").
        mutation_type:    One of: loss_of_function, gain_of_function, missense,
                          dominant_negative, haploinsufficiency, frameshift,
                          nonsense, deletion, amplification, general.
        specific_variant: Optional HGVS notation (e.g. "p.R175H", "c.185delAG").
        organism:         Host organism for taxid resolution (default Homo sapiens).
        step_cb:          Optional async progress callback.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    report = ReverseGeneticsReport(
        gene=gene.upper(),
        organism=organism,
        mutation_type=mutation_type,
        specific_variant=specific_variant,
    )

    async with httpx.AsyncClient(timeout=30) as client:

        # ── 1. UniProt — protein context ────────────────────────────────────
        await _step(f"[1/9] Fetching {gene} protein data from UniProt…")
        uniprot_entry: Optional[dict] = None
        ncbi_gene_id: Optional[str] = None
        try:
            from fetchers import uniprot as uni
            results = await uni.search_protein(
                f"reviewed:true AND gene_exact:{gene}",
                organism if organism != "Homo sapiens" else "Homo sapiens",
                size=1,
            )
            if results:
                acc = results[0]["primaryAccession"]
                uniprot_entry = await uni.fetch_by_accession(acc)
                report.uniprot_id = acc

            if isinstance(uniprot_entry, dict):
                report.protein_name = (
                    uniprot_entry.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "")
                ) or gene
                report.protein_function = _extract_protein_function(uniprot_entry)

                try:
                    from fetchers.brenda import get_ec_number
                    report.ec_number = get_ec_number(uniprot_entry)
                except Exception:
                    pass

                ncbi_gene_id = uni.get_ncbi_gene_id(uniprot_entry)
                report.is_essential, report.is_tumor_suppressor, report.is_oncogene = (
                    _infer_gene_role(uniprot_entry)
                )

                # Domain / feature impact analysis
                report.domain_impacts, report.target_domain = _analyze_domain_impacts(
                    uniprot_entry, specific_variant, mutation_type
                )
        except Exception:
            pass

        # ── 2. ClinVar — known pathogenic variants ──────────────────────────
        await _step(f"[2/9] Querying ClinVar for {gene} variants…")
        try:
            from fetchers.ncbi_clinvar import fetch_all_clinvar_entries
            clinvar_entries = await fetch_all_clinvar_entries(gene)
            report.known_clinvar_variants = clinvar_entries[:20]
            for e in clinvar_entries:
                sig = e.get("significance", "").lower()
                if "pathogenic" in sig or "likely pathogenic" in sig:
                    report.clinical_correlations.append(ClinicalCorrelation(
                        disease=e.get("disease", "Unknown"),
                        significance=e.get("significance", "Unknown"),
                        source="ClinVar",
                    ))
        except Exception:
            pass

        # ── 3. OMIM — mendelian disease associations ─────────────────────────
        await _step(f"[3/9] Fetching OMIM disease associations for {gene}…")
        try:
            from fetchers.omim import fetch_omim_diseases
            diseases, mim_ids = await fetch_omim_diseases(gene)
            report.omim_diseases = diseases[:8]
            for disease, mim_id in zip(diseases[:5], mim_ids[:5]):
                # Avoid duplicating diseases already in ClinVar correlations
                if not any(c.disease == disease for c in report.clinical_correlations):
                    report.clinical_correlations.append(ClinicalCorrelation(
                        disease=disease,
                        significance="Pathogenic (OMIM)",
                        omim_id=mim_id,
                        source="OMIM",
                    ))
        except Exception:
            pass

        # ── 4. KEGG — metabolic pathway disruption ──────────────────────────
        await _step(f"[4/9] Mapping disrupted KEGG metabolic pathways…")
        if ncbi_gene_id:
            try:
                report.disrupted_kegg_pathways = await _fetch_kegg_pathways_for_gene(
                    ncbi_gene_id, client
                )
                report.pathway_consequences = await _build_pathway_consequences(
                    report.disrupted_kegg_pathways, mutation_type, gene
                )
            except Exception:
                pass

        # ── 5. Reactome — biological process disruption ─────────────────────
        await _step(f"[5/9] Fetching disrupted Reactome pathways…")
        if report.uniprot_id:
            try:
                from fetchers.reactome import fetch_pathways
                reactome_pws = await fetch_pathways(report.uniprot_id)
                report.disrupted_reactome_pathways = [
                    {"id": p.get("stId", ""), "name": p.get("displayName", ""),
                     "species": p.get("speciesName", "")}
                    for p in reactome_pws[:10]
                ]
                # Add Reactome consequences
                for pw in reactome_pws[:5]:
                    report.pathway_consequences.append(PathwayConsequence(
                        pathway_id=pw.get("stId", ""),
                        pathway_name=pw.get("displayName", "?"),
                        consequence_type="disrupted" if "loss" in mutation_type.lower() else "altered",
                        downstream_effects=[
                            f"{gene} participates in this pathway — mutation alters its contribution",
                            "Dependent downstream reactions may be impaired",
                        ],
                        severity="Moderate",
                        source="Reactome",
                    ))
            except Exception:
                pass

        # ── 6. STRING — PPI network impact ──────────────────────────────────
        await _step(f"[6/9] Analysing PPI network disruption (STRING)…")
        try:
            ppi_partners = await _fetch_ppi_partners(gene, _HUMAN_TAXID, client)
            report.lost_interactions = ppi_partners[:15]
            hub_label, hub_score = _classify_hub(len(ppi_partners))
            report.hub_classification = hub_label
            report.hub_score = hub_score
        except Exception:
            pass

        # ── 7. PharmGKB — drug response alterations ──────────────────────────
        await _step(f"[7/9] Checking PharmGKB drug response implications…")
        try:
            from fetchers.pharmgkb import get_pharmgkb_id, fetch_drug_interactions, fetch_pgx_variants
            pgkb_id = await get_pharmgkb_id(gene)
            if pgkb_id:
                drug_ints, pgx_vars = await asyncio.gather(
                    fetch_drug_interactions(pgkb_id),
                    fetch_pgx_variants(pgkb_id),
                    return_exceptions=True,
                )
                if isinstance(drug_ints, list):
                    report.pgx_drug_interactions = drug_ints[:10]
                if isinstance(pgx_vars, list):
                    report.pgx_variants = pgx_vars[:10]
                    for v in pgx_vars[:5]:
                        report.clinical_correlations.append(ClinicalCorrelation(
                            disease=f"Drug response: {v.get('drug', '?')}",
                            significance=v.get("significance", "PGx variant"),
                            source="PharmGKB",
                        ))
        except Exception:
            pass

        # ── 8. Paralog discovery (compensatory potential) ────────────────────
        await _step(f"[8/9] Searching for compensatory paralogs…")
        try:
            report.paralogs = await _find_paralogs(gene, client)
        except Exception:
            pass

    # ── 9. Gemini synthesis ──────────────────────────────────────────────────
    await _step(f"[9/9] Gemini 2.5 Pro synthesizing reverse genetics analysis…")
    report = await _synthesize_with_gemini(report)

    report.references = _RG_REFERENCES
    return report
