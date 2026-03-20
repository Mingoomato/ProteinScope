"""Query orchestrator (query_engine.py).

Runs all fetchers concurrently via asyncio.gather, passes results through
analyzers, assembles a ProteinRecord, and dispatches to the output writer.
"""

from __future__ import annotations

import asyncio
import sys
from pathlib import Path
from typing import Optional

# Ensure package root is on sys.path when run directly
_ROOT = str(Path(__file__).parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from core.cache import Cache
from core.models import (
    ClinicalEntry,
    CrossSpeciesEntry,
    Isoform,
    ProteinRecord,
    Reference,
    SequenceFeature,
)

from fetchers import uniprot as uni


class ProteinSuggestionsError(Exception):
    """Raised when a query can't resolve to one protein but suggestions are available."""
    def __init__(self, query: str, suggestions: list[dict], organism: str | None = None):
        self.query = query
        self.suggestions = suggestions
        self.organism = organism
        super().__init__(f"Query '{query}' matched multiple proteins.")


class ProteinCompatibilityRequest(Exception):
    """Raised when a query is detected as a cross-species substitution question."""
    def __init__(self, source_gene: str, source_organism: str,
                 target_organism: str, target_cell: str | None):
        self.source_gene = source_gene
        self.source_organism = source_organism
        self.target_organism = target_organism
        self.target_cell = target_cell
        super().__init__(f"Compatibility query: {source_gene} from {source_organism} in {target_organism}")


class TaxonAmbiguityError(Exception):
    """Raised when an organism name refers to a broad taxonomic group."""
    def __init__(self, organism: str, taxon_rank: str, species: list[dict]):
        self.organism = organism
        self.taxon_rank = taxon_rank
        self.species = species
        super().__init__(f"'{organism}' is a broad {taxon_rank}-level taxon; species selection required.")


class ReverseGeneticsRequest(Exception):
    """Raised when a query asks about the effect of a gene mutation (reverse genetics)."""
    def __init__(self, gene: str, mutation_type: str,
                 specific_variant: Optional[str] = None, organism: str = "Homo sapiens"):
        self.gene = gene
        self.mutation_type = mutation_type
        self.specific_variant = specific_variant
        self.organism = organism
        super().__init__(f"Reverse genetics: {gene} ({mutation_type})")


class RNAiRequest(Exception):
    """Raised when a query asks to design or analyze an RNAi experiment."""
    def __init__(self, gene: str, rnai_type: str = "siRNA",
                 cell_type: Optional[str] = None, organism: str = "Homo sapiens"):
        self.gene = gene
        self.rnai_type = rnai_type
        self.cell_type = cell_type
        self.organism = organism
        super().__init__(f"RNAi: {rnai_type} knockdown of {gene}")


class PathwayAnalysisRequest(Exception):
    """Raised when a query asks to engineer a complete multi-protein biological pathway."""
    def __init__(self, query: str, plan: dict):
        self.query = query
        self.plan = plan
        super().__init__(f"Pathway engineering: {plan.get('pathway_name', query)}")


class OffTopicQueryError(Exception):
    """Raised when a query has no biological/protein context ProteinScope can process."""
    def __init__(self, query: str, reason: str = ""):
        self.query = query
        self.reason = reason
        super().__init__(f"Off-topic query: '{query}' — {reason}")


class SpeciesFallbackError(Exception):
    """Raised when UniProt found no results for the requested organism and fell back to another species.

    Carries the fallback accession so the caller can offer the user a choice to proceed or cancel.
    """
    def __init__(self, accession: str, requested_organism: str, actual_organism: str, gene: str = ""):
        self.accession = accession
        self.requested_organism = requested_organism
        self.actual_organism = actual_organism
        self.gene = gene
        super().__init__(
            f"No '{gene or accession}' found in {requested_organism}; "
            f"fell back to {actual_organism} ({accession})."
        )


from fetchers.ncbi_gene import fetch_cds_for_gene
from fetchers.ncbi_clinvar import fetch_all_clinvar_entries
from fetchers.omim import fetch_omim_diseases
from fetchers.brenda import (
    fetch_brenda_temperature,
    fetch_brenda_ph,
    extract_temperature_from_uniprot,
    extract_ph_from_uniprot,
    get_ec_number,
)
from fetchers.alphafold import fetch_alphafold
from fetchers.pdb import fetch_experimental_structures
from fetchers.pubmed import fetch_citations

from analyzers.domain_extractor import extract_features_from_uniprot
from analyzers.isoform_handler import fetch_isoform_fastas, build_isoforms
from analyzers.cross_species import analyze_cross_species
from analyzers.reference_builder import populate_citations
from analyzers.pgx_reporter import build_drug_interactions, build_pgx_variants
from analyzers.pathway_integrator import (
    build_metabolic_pathways_from_kegg,
    build_signaling_pathways_from_reactome,
    build_disease_pathways_from_reactome,
    build_pharmacogenomic_pathways_from_pharmgkb,
    build_protein_interactions_from_string,
)


def _build_references(raw_refs: list[dict]) -> list[Reference]:
    refs = []
    for r in raw_refs:
        try:
            year_val = int(str(r.get("year", "0"))[:4])
        except ValueError:
            year_val = 0
        refs.append(Reference(
            pubmed_id=r.get("pubmed_id"),
            doi=r.get("doi"),
            title=r.get("title", ""),
            authors=r.get("authors", []),
            journal=r.get("journal", ""),
            year=year_val,
        ))
    return refs


async def _resolve_accession(query: str, organism: Optional[str], session_context: str = "") -> str:
    """Return a UniProt accession for the query string.

    If the query looks like a UniProt accession (e.g. P68871), use it directly.
    Otherwise search UniProt and return the top hit's accession.
    """
    # UniProt accession format (from UniProt docs):
    #   [OPQ][0-9][A-Z0-9]{3}[0-9]                    e.g. P68871, Q9Y6K9
    #   [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}    e.g. A0A000, B2RUZ4
    # Position 2 (0-indexed) is always a digit — this excludes gene symbols like CYP2C9.
    import re
    if re.match(
        r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$",
        query.upper(),
    ):
        return query.upper()

    # ---------------------------------------------------------------------------
    # Step 1 — NL pipeline (always runs first to understand query intent)
    # ---------------------------------------------------------------------------
    # Gemini classifies the query BEFORE any UniProt search so that:
    #   • Multi-protein / pathway queries → ProteinSuggestionsError (suggestion cards)
    #   • Cross-species queries           → ProteinCompatibilityRequest (compat flow)
    #   • Broad taxon queries             → TaxonAmbiguityError (species selection)
    #   • Single-protein queries          → extracts gene symbol + organism, then searches
    #
    # Simple gene symbols (e.g. "POLD1") are detected and fast-pathed through direct
    # UniProt search AFTER NL, not before, so intent is always checked first.
    # ---------------------------------------------------------------------------

    try:
        from core.gemini_interpreter import (
            classify_query_intent,
            parse_nl_protein_query,
            resolve_nl_query,
            suggest_proteins_for_query,
            translate_to_english,
        )

        english_query = await translate_to_english(query)

        # --- Step 0: Top-level intent classification ---
        # Runs before anything else so that off-topic / vague meta-commands
        # (e.g. "Make a report suggesting why this is impossible.") are caught
        # immediately rather than producing nonsensical protein suggestions.
        intent_result = await classify_query_intent(english_query, session_context=session_context)
        if intent_result["intent"] == "off_topic":
            raise OffTopicQueryError(query, intent_result.get("reason", ""))

        # --- Pathway engineering check (before single-protein compatibility) ---
        # Catches "enable photosynthesis in human cells" BEFORE it falls through to
        # the compatibility check which would only pick up one protein (e.g. RBCL).
        try:
            from core.gemini_interpreter import parse_pathway_analysis_plan
            pathway = await parse_pathway_analysis_plan(english_query)
            if pathway["is_pathway_analysis"]:
                raise PathwayAnalysisRequest(query, pathway)
        except PathwayAnalysisRequest:
            raise
        except Exception:
            pass

        # --- Compatibility check ---
        try:
            from core.gemini_interpreter import parse_compatibility_query
            compat = await parse_compatibility_query(english_query)
            if compat["is_compatibility"] and compat["source_organism"]:
                source_gene = compat["source_gene"] or await resolve_nl_query(english_query)
                if source_gene:
                    raise ProteinCompatibilityRequest(
                        source_gene=source_gene,
                        source_organism=compat["source_organism"],
                        target_organism=compat["target_organism"] or "Homo sapiens",
                        target_cell=compat["target_cell"],
                    )
        except ProteinCompatibilityRequest:
            raise
        except Exception:
            pass

        # --- RNAi experiment design check (before reverse genetics — more specific) ---
        try:
            from core.gemini_interpreter import parse_rnai_query
            rnai = await parse_rnai_query(english_query)
            if rnai["is_rnai"] and rnai["gene"]:
                raise RNAiRequest(
                    gene=rnai["gene"],
                    rnai_type=rnai.get("rnai_type", "siRNA"),
                    cell_type=rnai.get("cell_type"),
                    organism=rnai.get("organism") or organism or "Homo sapiens",
                )
        except RNAiRequest:
            raise
        except Exception:
            pass

        # --- Reverse genetics / mutation impact check ---
        try:
            from core.gemini_interpreter import parse_reverse_genetics_query
            rg = await parse_reverse_genetics_query(english_query)
            if rg["is_reverse_genetics"] and rg["gene"]:
                raise ReverseGeneticsRequest(
                    gene=rg["gene"],
                    mutation_type=rg.get("mutation_type", "general"),
                    specific_variant=rg.get("specific_variant"),
                    organism=rg.get("organism") or organism or "Homo sapiens",
                )
        except ReverseGeneticsRequest:
            raise
        except Exception:
            pass

        # --- Intent + organism extraction ---
        # Always parse the ORIGINAL query first — parse_nl_protein_query explicitly
        # supports any language (Korean, Japanese, etc.), so "팽귄" → "penguin" is
        # correctly understood without going through the English translation that
        # may silently drop the organism name.
        parsed = await parse_nl_protein_query(query, session_context=session_context)
        gemini_organism = parsed.get("organism")

        # If the original-query parse didn't yield an organism, try the English
        # translation as a second attempt and merge any missing fields.
        if not gemini_organism and query != english_query:
            try:
                parsed_en = await parse_nl_protein_query(english_query, session_context=session_context)
                gemini_organism = parsed_en.get("organism")
                if parsed_en.get("is_multi"):
                    parsed["is_multi"] = True
                if not parsed.get("gene") and parsed_en.get("gene"):
                    parsed["gene"] = parsed_en["gene"]
            except Exception:
                pass

        # Anti-hallucination check — validate against the ENGLISH translation.
        # If an organism is legitimately mentioned (in any language), it will appear
        # in the English translation: "팽귄" / "פינגווין" → "penguin" (in english_query).
        # A hallucinated "human" will NOT appear in "DNA replication proteins in penguin".
        # This catches hallucinations across all input languages uniformly.
        if gemini_organism and gemini_organism.lower() not in english_query.lower():
            gemini_organism = None

        extracted_organism = organism or gemini_organism  # user-supplied always wins

        # --- Taxon ambiguity check (only when user hasn't specified a precise species) ---
        if extracted_organism and not organism:
            try:
                from core.gemini_interpreter import detect_taxon_ambiguity
                taxon = await detect_taxon_ambiguity(extracted_organism)
                if taxon["is_broad"] and taxon["species"]:
                    raise TaxonAmbiguityError(
                        organism=extracted_organism,
                        taxon_rank=taxon["taxon_rank"],
                        species=taxon["species"],
                    )
            except TaxonAmbiguityError:
                raise
            except Exception:
                pass

        # --- Multi-protein keyword heuristic override ---
        _MULTI_KEYWORDS = (
            " proteins", " enzymes", " factors", " components", " all ", " list ",
            "pathway", "process", "complex", "involved in", "related to",
            "associated with", "participating in",
        )
        if not parsed["is_multi"] and any(kw in english_query.lower() for kw in _MULTI_KEYWORDS):
            parsed["is_multi"] = True

        # --- Meta-command guard: block "make a report / why / explain" with no protein context ---
        # These are instructions, not protein queries. Without a gene or organism they would
        # produce meaningless protein suggestions (e.g. "Make a report" → apoptosis proteins).
        _META_PREFIXES = (
            "make ", "create ", "generate ", "write ", "explain ", "describe ",
            "why ", "how ", "tell me", "show me", "give me",
        )
        _eq_lower = english_query.lower().strip()
        _is_meta = (
            any(_eq_lower.startswith(p) for p in _META_PREFIXES)
            and not parsed.get("gene")
            and not extracted_organism
        )
        if _is_meta:
            parsed["is_multi"] = False  # don't trigger protein suggestions for meta-commands

        # --- Route: multi-protein vs single-protein ---
        if parsed["is_multi"]:
            suggestions = await suggest_proteins_for_query(english_query, organism=extracted_organism)
            if suggestions:
                raise ProteinSuggestionsError(query, suggestions, organism=extracted_organism)

        else:
            gene = parsed.get("gene") or await resolve_nl_query(english_query)
            if gene:
                found_org = None
                found_results = None
                for org_attempt in ([extracted_organism, None] if extracted_organism else [None]):
                    nl_results = await uni.search_protein(gene, org_attempt, size=1)
                    if nl_results:
                        found_org = org_attempt
                        found_results = nl_results
                        break
                if found_results:
                    accession = found_results[0]["primaryAccession"]
                    if extracted_organism and found_org is None:
                        from fetchers.uniprot import extract_organism as _xorg, extract_gene_name as _xgene
                        actual_org, _ = _xorg(found_results[0])
                        raise SpeciesFallbackError(
                            accession=accession,
                            requested_organism=extracted_organism,
                            actual_organism=actual_org or "unknown species",
                            gene=gene,
                        )
                    return accession

    except (ProteinSuggestionsError, ProteinCompatibilityRequest, PathwayAnalysisRequest,
            TaxonAmbiguityError, SpeciesFallbackError, OffTopicQueryError,
            ReverseGeneticsRequest, RNAiRequest):
        raise
    except Exception:
        pass

    # ---------------------------------------------------------------------------
    # Step 2 — Direct UniProt search (fast-path for plain gene symbols / names)
    # ---------------------------------------------------------------------------
    # Only reached when the NL pipeline couldn't resolve the query (e.g. the
    # query was already a bare gene symbol like "POLD1" or "CYP2C9" and Gemini
    # classified it as a single protein but UniProt returned nothing via NL path).
    # ---------------------------------------------------------------------------

    _is_simple_query = " " not in query.strip() and all(ord(c) < 128 for c in query)
    candidates = [
        f"reviewed:true AND gene_exact:{query}",
        f"reviewed:true AND protein_name:{query}",
        f"reviewed:true AND ({query})",
        query,
    ] if _is_simple_query else [
        f"reviewed:true AND gene_exact:{query}",
        f"reviewed:true AND protein_name:{query}",
    ]

    organism_attempts = [organism, None] if organism else [None]
    results = []
    found_with_org = None
    for org_try in organism_attempts:
        for q in candidates:
            results = await uni.search_protein(q, org_try, size=1)
            if results:
                found_with_org = org_try
                break
        if results:
            break

    if results:
        accession = results[0]["primaryAccession"]
        if organism and found_with_org is None:
            from fetchers.uniprot import extract_organism as _xorg, extract_gene_name as _xgene
            actual_org, _ = _xorg(results[0])
            raise SpeciesFallbackError(
                accession=accession,
                requested_organism=organism,
                actual_organism=actual_org or "unknown species",
                gene=_xgene(results[0]) or query,
            )
        return accession

    raise ValueError(f"No UniProt results found for query: {query!r}")


async def fetch_all(
    accession: str,
    organism: Optional[str],
    cross_species: bool,
    cache: Optional[Cache],
    include_pathways: bool = True,
    no_images: bool = False,
    min_pgx_evidence: str = "4",
    interaction_score: float = 0.4,
    ai_summary: bool = False,
    no_annotate: bool = False,
    max_callouts: int = 10,
    no_ai_polish: bool = False,
    session_context: str = "",
) -> ProteinRecord:
    """Run all fetchers concurrently and assemble a ProteinRecord."""

    # Ph4: resolve session entity aliases before lookup
    if session_context:
        try:
            from core.session_manager import session_manager as _sm
            # session_context is a pre-built block; entity map isn't directly accessible here
            # The resolution was done in app.py before calling fetch_all()
            pass
        except Exception:
            pass

    # ── UniProt (always needed first to get IDs for other fetchers) ──
    cached_entry = cache.get(accession, "uniprot") if cache else None
    if cached_entry:
        entry = cached_entry
    else:
        entry = await uni.fetch_by_accession(accession)
        if cache:
            cache.set(accession, "uniprot", entry)

    # ── Extract identity fields ──
    protein_name = uni.extract_protein_name(entry)
    gene_name = uni.extract_gene_name(entry)
    organism_name, taxonomy_id = uni.extract_organism(entry)
    canonical_sequence = entry.get("sequence", {}).get("value", "")
    sequence_length = entry.get("sequence", {}).get("length", len(canonical_sequence))
    ncbi_gene_id = uni.get_ncbi_gene_id(entry)

    # ── Extract features from UniProt ──
    features_raw = entry.get("features", [])
    binding_sites, active_sites, domains = extract_features_from_uniprot(
        features_raw, canonical_sequence
    )

    # ── Gather remaining data concurrently ──
    ec_number = get_ec_number(entry)
    effective_organism = organism or organism_name

    async def _noop_cds():
        return (None, False)

    async def _noop():
        return None

    tasks = [
        fetch_cds_for_gene(ncbi_gene_id) if ncbi_gene_id else _noop_cds(),
        fetch_all_clinvar_entries(gene_name),
        fetch_omim_diseases(gene_name),
        fetch_brenda_temperature(ec_number, effective_organism) if ec_number else _noop(),
        fetch_brenda_ph(ec_number, effective_organism) if ec_number else _noop(),
        fetch_alphafold(accession),
        fetch_experimental_structures(accession),
        fetch_isoform_fastas(accession),
    ]

    (
        cds_result,
        clinvar_entries,
        omim_result,
        temp_brenda,
        ph_brenda,
        alphafold_result,
        pdb_structures,
        isoform_seqs,
    ) = await asyncio.gather(*tasks, return_exceptions=True)

    # ── Unpack CDS ──
    if isinstance(cds_result, tuple):
        cds_seq, is_spliced = cds_result
    else:
        cds_seq, is_spliced = None, False

    # ── Unpack OMIM ──
    if isinstance(omim_result, tuple):
        omim_diseases, omim_deficiencies = omim_result
    else:
        omim_diseases, omim_deficiencies = [], []

    # ── Unpack AlphaFold ──
    if isinstance(alphafold_result, tuple):
        af_url, af_plddt = alphafold_result
    else:
        af_url, af_plddt = None, None

    # ── PDB IDs ──
    pdb_ids = []
    if isinstance(pdb_structures, list):
        pdb_ids = [s["pdb_id"] for s in pdb_structures]

    # ── Temperature / pH ──
    temp = temp_brenda if isinstance(temp_brenda, float) else extract_temperature_from_uniprot(entry)
    ph = ph_brenda if isinstance(ph_brenda, float) else extract_ph_from_uniprot(entry)

    # ── Isoforms ──
    isoforms: list[Isoform] = []
    if isinstance(isoform_seqs, dict):
        isoforms = build_isoforms(entry, canonical_sequence, isoform_seqs)

    # ── Clinical variants ──
    clinical_variants: list[ClinicalEntry] = []
    if isinstance(clinvar_entries, list):
        for cv in clinvar_entries:
            from core.evidence import DataProvenance, EvidenceGrade
            clinical_variants.append(ClinicalEntry(
                variant=cv.get("variant"),
                disease=cv.get("disease", ""),
                significance=cv.get("significance", "unknown"),
                provenance=DataProvenance(
                    source="ClinVar (NCBI REST v2)",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                ),
            ))

    # ── P1: ESM-2 variant fitness scoring ──
    # Parse HGVS p.A123V notation from clinical_variants and score them.
    # Runs synchronously in executor to avoid blocking the event loop.
    variant_fitness_scores = []
    if canonical_sequence and clinical_variants:
        try:
            import re as _re
            from analyzers.variant_scorer import score_variants
            _hgvs_pat = _re.compile(r"p\.([A-Z])(\d+)([A-Z])", _re.IGNORECASE)
            to_score: list[tuple[str, int, str]] = []
            for cv in clinical_variants:
                v_str = cv.variant or ""
                m = _hgvs_pat.search(v_str)
                if m:
                    ref_aa, pos_str, alt_aa = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
                    if 1 <= pos_str <= len(canonical_sequence):
                        to_score.append((ref_aa, pos_str, alt_aa))
            if to_score:
                loop = asyncio.get_event_loop()
                scored = await loop.run_in_executor(
                    None, score_variants, canonical_sequence, to_score
                )
                # Annotate ClinicalEntry objects in-place
                score_map = {f"{s.ref_aa}{s.position}{s.alt_aa}": s for s in scored}
                for cv in clinical_variants:
                    v_str = cv.variant or ""
                    m = _hgvs_pat.search(v_str)
                    if m:
                        key = f"{m.group(1).upper()}{int(m.group(2))}{m.group(3).upper()}"
                        s = score_map.get(key)
                        if s and s.delta_log_likelihood is not None:
                            cv.esm2_score = s.delta_log_likelihood
                            cv.esm2_interpretation = s.predicted_impact
                variant_fitness_scores = scored
        except Exception:
            pass  # fair-esm not installed or model load failed — degrade gracefully

    # ── P2: PTM Logic Engine ──
    ptm_logic_result = None
    try:
        from fetchers.phosphosite import fetch_ptm_sites
        from analyzers.ptm_logic_engine import run_ptm_analysis
        ptm_raw = await fetch_ptm_sites(gene_name, entry)
        if ptm_raw:
            ptm_logic_result = await run_ptm_analysis(gene_name, ptm_raw)
    except Exception:
        pass

    # ── P3: IDP / LLPS Detection (auto-triggers when pLDDT < 70) ──
    idp_analysis_result = None
    if canonical_sequence and af_plddt is not None and af_plddt < 70:
        try:
            from fetchers.iupred3 import fetch_iupred_scores
            from fetchers.fuzdrop import fetch_fuzdrop
            from analyzers.idp_analyzer import run_idp_analysis
            iupred_data, fuzdrop_data = await asyncio.gather(
                fetch_iupred_scores(canonical_sequence),
                fetch_fuzdrop(canonical_sequence),
                return_exceptions=True,
            )
            iupred_data = iupred_data if isinstance(iupred_data, dict) else {}
            fuzdrop_data = fuzdrop_data if isinstance(fuzdrop_data, dict) else {}
            idp_analysis_result = await run_idp_analysis(
                gene=gene_name,
                sequence=canonical_sequence,
                iupred_data=iupred_data,
                fuzdrop_data=fuzdrop_data,
                af_plddt=af_plddt,
            )
        except Exception:
            pass

    # ── P4: Metalloenzyme + PROPKA (auto-triggers when metal cofactors detected) ──
    metalloenzyme_result = None
    _METALS = {"zinc", "iron", "copper", "manganese", "magnesium", "calcium", "cobalt",
               "nickel", "molybdenum", "vanadium", "tungsten"}
    _cofactors_raw = uni.extract_cofactors(entry)
    _cofactor_str = " ".join(_cofactors_raw).lower()
    if any(m in _cofactor_str for m in _METALS):
        try:
            from analyzers.metalloenzyme_analyzer import run_metalloenzyme_analysis
            metalloenzyme_result = await run_metalloenzyme_analysis(
                gene=gene_name,
                cofactors=_cofactors_raw,
                uniprot_entry=entry,
                af_pdb_url=af_url,
            )
        except Exception:
            pass

    # ── P5: Host-Protein Expression Compatibility Matrix (pure computation, no I/O) ──
    expression_compatibility_result = None
    if canonical_sequence:
        try:
            from analyzers.host_compatibility import run_expression_compatibility
            expression_compatibility_result = run_expression_compatibility(
                gene=gene_name,
                canonical_sequence=canonical_sequence,
                encoding_dna_cds=cds_seq,
                entry_features=entry.get("features", []),
            )
        except Exception:
            pass

    # ── P6: Binding Pocket Druggability (conditional on AlphaFold URL) ──
    druggability_result = None
    if af_url:
        try:
            from analyzers.druggability_analyzer import run_druggability_analysis
            druggability_result = await run_druggability_analysis(
                gene=gene_name,
                af_pdb_url=af_url,
                uniprot_entry=entry,
            )
        except Exception:
            pass

    # ── Ph2: Conformational ANM/GNM (conditional on AlphaFold URL) ──
    conformational_result = None
    if af_url:
        try:
            from analyzers.conformational_analyzer import run_conformational_analysis
            conformational_result = await run_conformational_analysis(
                gene=gene_name,
                af_pdb_url=af_url,
            )
        except Exception:
            pass

    # ── P9: Inverse Folding via ESM-IF1 (conditional on AlphaFold URL) ──
    inverse_designs_result = []
    if af_url:
        try:
            from analyzers.inverse_folder import design_sequences
            inverse_designs_result = await design_sequences(
                gene=gene_name,
                af_pdb_url=af_url,
                af_plddt=af_plddt,
            )
        except Exception:
            pass

    # ── P10: Therapeutic Decision Simulator ──
    therapeutic_result = None
    try:
        from fetchers.clinicaltrials import fetch_clinical_trials
        from fetchers.protein_atlas import fetch_hpa_protein_expression
        from analyzers.therapeutic_simulator import run_therapeutic_analysis
        _ct_data, _hpa_data = await asyncio.gather(
            fetch_clinical_trials(gene_name),
            fetch_hpa_protein_expression(gene_name),
            return_exceptions=True,
        )
        _ct_data = _ct_data if isinstance(_ct_data, list) else []
        _hpa_data = _hpa_data if isinstance(_hpa_data, dict) else {}
        therapeutic_result = await run_therapeutic_analysis(
            gene=gene_name,
            subcellular_locations=uni.extract_subcellular_location(entry),
            function_description=uni.extract_function(entry),
            clinical_trials=_ct_data,
            protein_atlas_data=_hpa_data if _hpa_data else None,
        )
    except Exception:
        pass

    # ── References ──
    raw_refs = uni.extract_references(entry)
    pmids = [r["pubmed_id"] for r in raw_refs if r.get("pubmed_id")]
    pubmed_data = await fetch_citations(pmids)

    # Merge PubMed metadata into raw_refs
    pmid_map = {r["pubmed_id"]: r for r in pubmed_data}
    merged_refs: list[dict] = []
    for ref in raw_refs:
        pm = pmid_map.get(ref.get("pubmed_id", ""), {})
        merged_refs.append({
            "pubmed_id": ref.get("pubmed_id"),
            "doi": pm.get("doi") or ref.get("doi"),
            "title": pm.get("title") or ref.get("title", ""),
            "authors": pm.get("authors") or ref.get("authors", []),
            "journal": pm.get("journal") or ref.get("journal", ""),
            "year": pm.get("year") or ref.get("year", 0),
        })
    references = populate_citations(_build_references(merged_refs))

    # ── Pathway data (concurrent with nothing else left to wait for) ──
    drug_interactions_list = []
    pgx_variants_list = []
    pharmacogenomic_pathways_list = []
    metabolic_pathways_list = []
    signaling_pathways_list = []
    protein_interactions_list = []
    disease_pathways_list = []

    if include_pathways:
        from fetchers.pharmgkb import fetch_all_pharmgkb
        from fetchers.dgidb import fetch_dgidb_interactions
        from fetchers.chembl import fetch_all_chembl
        from fetchers.kegg import fetch_all_kegg
        from fetchers.metacyc import fetch_all_metacyc
        from fetchers.reactome import fetch_all_reactome
        from fetchers.string_db import fetch_all_string

        (
            pharmgkb_data,
            dgidb_data,
            chembl_data,
            kegg_data,
            _metacyc_data,
            reactome_data,
            string_data,
        ) = await asyncio.gather(
            fetch_all_pharmgkb(gene_name),
            fetch_dgidb_interactions(gene_name),
            fetch_all_chembl(accession),
            fetch_all_kegg(ncbi_gene_id or ""),
            fetch_all_metacyc(gene_name),
            fetch_all_reactome(accession),
            fetch_all_string(gene_name, min_score=int(interaction_score * 1000)),
            return_exceptions=True,
        )

        if isinstance(pharmgkb_data, Exception):
            pharmgkb_data = {}
        if isinstance(dgidb_data, Exception):
            dgidb_data = []
        if isinstance(chembl_data, Exception):
            chembl_data = {}
        if isinstance(kegg_data, Exception):
            kegg_data = {}
        if isinstance(reactome_data, Exception):
            reactome_data = {}
        if isinstance(string_data, Exception):
            string_data = {}

        drug_interactions_list = build_drug_interactions(
            pharmgkb_data, dgidb_data if isinstance(dgidb_data, list) else [],
            chembl_data, min_pgx_evidence=min_pgx_evidence,
        )
        pgx_variants_list = build_pgx_variants(pharmgkb_data, gene_name, min_pgx_evidence)
        pharmacogenomic_pathways_list = build_pharmacogenomic_pathways_from_pharmgkb(pharmgkb_data)
        metabolic_pathways_list = build_metabolic_pathways_from_kegg(kegg_data, gene_name)
        signaling_pathways_list = build_signaling_pathways_from_reactome(reactome_data, accession, gene_name)
        disease_pathways_list = build_disease_pathways_from_reactome(reactome_data, kegg_data, gene_name)
        protein_interactions_list = build_protein_interactions_from_string(string_data, min_score=interaction_score)

    # ── Layer 1: Allosteric Network ──
    allosteric_result = None
    if af_url:
        try:
            from fetchers.phasepdb import fetch_phasepdb
            from analyzers.allosteric_analyzer import run_allosteric_analysis
            _phasepdb_data = await fetch_phasepdb(accession)
            allosteric_result = await run_allosteric_analysis(
                gene=gene_name, af_pdb_url=af_url, phasepdb_data=_phasepdb_data,
            )
        except Exception:
            pass

    # ── Layer 1: Epistasis Analysis ──
    epistasis_result = None
    if canonical_sequence and variant_fitness_scores:
        try:
            from fetchers.mavedb import fetch_mavedb
            from analyzers.epistasis_analyzer import run_epistasis_analysis
            _mavedb_data = await fetch_mavedb(gene_name)
            epistasis_result = await run_epistasis_analysis(
                gene=gene_name, sequence=canonical_sequence,
                variant_fitness_scores=variant_fitness_scores, mavedb_data=_mavedb_data,
            )
        except Exception:
            pass

    # ── Layer 1: Binding Energy Estimation ──
    binding_energy_result = None
    if drug_interactions_list:
        try:
            from analyzers.binding_energy_estimator import run_binding_energy_estimation
            binding_energy_result = await run_binding_energy_estimation(
                gene=gene_name, drug_interactions=drug_interactions_list,
                druggability_report=druggability_result,
            )
        except Exception:
            pass

    # ── Layer 1: Observable Prediction (NMR/FRET/HDX) ──
    observables_result = None
    if canonical_sequence:
        try:
            from analyzers.observable_predictor import run_observable_prediction
            observables_result = await run_observable_prediction(
                gene=gene_name, sequence=canonical_sequence, af_plddt=af_plddt,
            )
        except Exception:
            pass

    # ── Layer 2: Active Learning Mutation Advisor ──
    active_learning_result = None
    if canonical_sequence and variant_fitness_scores:
        try:
            from fetchers.mavedb import fetch_mavedb as _fetch_mavedb2
            from analyzers.active_learning_advisor import run_active_learning_advice
            _mavedb2 = await _fetch_mavedb2(gene_name)
            active_learning_result = await run_active_learning_advice(
                gene=gene_name, sequence=canonical_sequence,
                variant_fitness_scores=variant_fitness_scores, mavedb_data=_mavedb2,
            )
        except Exception:
            pass

    # ── Layer 3: Glycan Analysis ──
    glycan_result = None
    if canonical_sequence:
        try:
            from fetchers.glycan import fetch_glycan_data
            from analyzers.glycan_analyzer import run_glycan_analysis
            _glycan_data = await fetch_glycan_data(accession)
            glycan_result = await run_glycan_analysis(
                gene=gene_name, sequence=canonical_sequence,
                entry_features=entry.get("features", []),
                glycan_data=_glycan_data,
            )
        except Exception:
            pass

    # ── Layer 3: ADMET Prediction ──
    admet_result = None
    if drug_interactions_list:
        try:
            from fetchers.admet import fetch_admet_profiles
            from analyzers.admet_analyzer import run_admet_analysis
            _admet_data = await fetch_admet_profiles(
                [{"name": getattr(d, "drug_name", ""), "drug_id": getattr(d, "drug_id", "")}
                 for d in drug_interactions_list]
            )
            admet_result = await run_admet_analysis(
                gene=gene_name, drug_interactions=drug_interactions_list,
                admet_data=_admet_data if isinstance(_admet_data, list) else [],
            )
        except Exception:
            pass

    # ── Layer 3: Directed Evolution Campaign ──
    directed_evolution_result = None
    if canonical_sequence:
        try:
            from analyzers.directed_evolution_analyzer import run_directed_evolution_design
            directed_evolution_result = await run_directed_evolution_design(
                gene=gene_name, sequence=canonical_sequence,
                variant_fitness_scores=variant_fitness_scores,
                binding_sites=binding_sites, active_sites=active_sites,
            )
        except Exception:
            pass

    # ── Layer 3: PROTAC Feasibility ──
    protac_result = None
    if canonical_sequence:
        try:
            from analyzers.protac_analyzer import run_protac_analysis
            protac_result = await run_protac_analysis(
                gene=gene_name, sequence=canonical_sequence,
                domains=domains, druggability_report=druggability_result,
            )
        except Exception:
            pass

    # ── Layer 3: Selectivity Landscape ──
    selectivity_result = None
    if druggability_result:
        try:
            from analyzers.selectivity_analyzer import run_selectivity_analysis
            selectivity_result = await run_selectivity_analysis(
                gene=gene_name, druggability_report=druggability_result,
                protein_interactions=protein_interactions_list,
            )
        except Exception:
            pass

    # ── Layer 3: Fragment Hotspot Map ──
    fragment_hotspot_result = None
    if af_url:
        try:
            from analyzers.fragment_hotspot_analyzer import run_fragment_hotspot_analysis
            fragment_hotspot_result = await run_fragment_hotspot_analysis(
                gene=gene_name, af_pdb_url=af_url,
                binding_sites=binding_sites, active_sites=active_sites,
            )
        except Exception:
            pass

    # ── Layer 3: AAV Gene Therapy Designer ──
    aav_result = None
    if canonical_sequence:
        try:
            from analyzers.aav_designer import run_aav_design
            aav_result = await run_aav_design(
                gene=gene_name, sequence=canonical_sequence,
                cds=cds_seq, function_description=uni.extract_function(entry),
            )
        except Exception:
            pass

    # ── Layer 3: CRISPR Integration Planner ──
    crispr_result = None
    if canonical_sequence:
        try:
            from analyzers.crispr_designer import run_crispr_design
            crispr_result = await run_crispr_design(
                gene=gene_name, sequence=canonical_sequence,
                cds=cds_seq, clinical_variants=clinical_variants,
            )
        except Exception:
            pass

    # ── Layer 3: Genetic Circuit Architect ──
    genetic_circuit_result = None
    if canonical_sequence:
        try:
            from analyzers.genetic_circuit_architect import run_genetic_circuit_design
            genetic_circuit_result = await run_genetic_circuit_design(
                gene=gene_name, sequence=canonical_sequence,
                function_description=uni.extract_function(entry),
            )
        except Exception:
            pass

    # ── Layer 3+: Covalent Inhibitor Designer ──
    covalent_result = None
    if canonical_sequence:
        try:
            from analyzers.covalent_inhibitor_designer import run_covalent_inhibitor_design
            covalent_result = await run_covalent_inhibitor_design(
                gene=gene_name, sequence=canonical_sequence,
                binding_sites=binding_sites, active_sites=active_sites,
            )
        except Exception:
            pass

    # ── Layer 3+: Protein Engineering Strategy ──
    engineering_result = None
    if canonical_sequence:
        try:
            from analyzers.engineering_strategy_recommender import run_engineering_strategy
            engineering_result = await run_engineering_strategy(
                gene=gene_name, sequence=canonical_sequence,
                af_plddt=af_plddt, binding_sites=binding_sites,
            )
        except Exception:
            pass

    # ── Layer 3+: H-Bond Network ──
    hbond_result = None
    if af_url:
        try:
            from analyzers.hbond_network_analyzer import run_hbond_analysis
            hbond_result = await run_hbond_analysis(
                gene=gene_name, af_pdb_url=af_url,
            )
        except Exception:
            pass

    # ── V3: Disease Association Layer (OpenTargets + DisGeNET + ClinGen + HPO) ──
    disease_association_result = None
    try:
        from fetchers.opentargets import fetch_opentargets, extract_ensembl_id_from_uniprot
        from fetchers.disgenet import fetch_disgenet
        from fetchers.clingen import fetch_clingen
        from fetchers.hpo import fetch_hpo_terms
        from analyzers.disease_association_analyzer import run_disease_association_analysis

        _ensembl_id = extract_ensembl_id_from_uniprot(entry)
        _ot_data, _disgenet_data, _clingen_data, _hpo_data = await asyncio.gather(
            fetch_opentargets(_ensembl_id or ""),
            fetch_disgenet(gene_name),
            fetch_clingen(gene_name),
            fetch_hpo_terms(ncbi_gene_id or ""),
            return_exceptions=True,
        )
        disease_association_result = await run_disease_association_analysis(
            gene=gene_name,
            ot_data=_ot_data if isinstance(_ot_data, list) else [],
            disgenet_data=_disgenet_data if isinstance(_disgenet_data, list) else [],
            clingen_data=_clingen_data if isinstance(_clingen_data, list) else [],
            hpo_data=_hpo_data if isinstance(_hpo_data, list) else [],
            step_cb=step_cb if "step_cb" in dir() else None,
        )
    except Exception:
        pass

    # ── V3: ETP Mapper — Beratan-Onuchic (triggers for redox cofactor proteins) ──
    etp_result = None
    _REDOX_KEYWORDS = {"heme", "haem", "cytochrome", "iron-sulfur", "fe-s", "fes",
                       "flavin", "fad", "fmn", "copper", "molybdenum", "4fe-4s", "2fe-2s"}
    _cofactors_str = " ".join(uni.extract_cofactors(entry)).lower()
    _has_redox = any(kw in _cofactors_str for kw in _REDOX_KEYWORDS)
    if af_url and _has_redox:
        try:
            from analyzers.etp_mapper import run_etp_analysis
            etp_result = await run_etp_analysis(
                gene=gene_name,
                af_pdb_url=af_url,
                cofactors=uni.extract_cofactors(entry),
                fragment_hotspot_map=fragment_hotspot_result,
                step_cb=step_cb if "step_cb" in dir() else None,
            )
        except Exception:
            pass

    # ── Cross-species ──
    cross_species_entries: list[CrossSpeciesEntry] = []
    if cross_species and binding_sites:
        first_binding_frag = binding_sites[0].sequence_fragment
        try:
            cross_species_entries = await analyze_cross_species(
                domain_fragment=first_binding_frag,
                reference_fragment=first_binding_frag,
            )
        except Exception:
            pass

    # ── GO terms ──
    bp_terms, mf_terms = uni.extract_go_terms(entry)

    # ── Assemble ──
    record = ProteinRecord(
        query=accession,
        uniprot_id=accession,
        protein_name=protein_name,
        gene_name=gene_name,
        organism=organism_name,
        taxonomy_id=taxonomy_id,
        function_description=uni.extract_function(entry),
        biological_process=bp_terms,
        molecular_function=mf_terms,
        cofactors=uni.extract_cofactors(entry),
        optimal_temperature_celsius=temp,
        optimal_ph=ph,
        subcellular_location=uni.extract_subcellular_location(entry),
        canonical_sequence=canonical_sequence,
        sequence_length=sequence_length,
        encoding_dna_cds=cds_seq,
        is_from_spliced_mrna=is_spliced,
        isoforms=isoforms,
        binding_sites=binding_sites,
        active_sites=active_sites,
        domains=domains,
        clinical_variants=clinical_variants,
        deficiency_diseases=omim_deficiencies,
        alphafold_pdb_url=af_url,
        alphafold_plddt_score=af_plddt,
        experimental_pdb_ids=pdb_ids,
        cross_species_comparison=cross_species_entries,
        references=references,
        drug_interactions=drug_interactions_list,
        pgx_variants=pgx_variants_list,
        pharmacogenomic_pathways=pharmacogenomic_pathways_list,
        metabolic_pathways=metabolic_pathways_list,
        signaling_pathways=signaling_pathways_list,
        protein_interactions=protein_interactions_list,
        disease_pathways=disease_pathways_list,
        variant_fitness_scores=variant_fitness_scores,
        ptm_logic=ptm_logic_result,
        idp_analysis=idp_analysis_result,
        metalloenzyme_analysis=metalloenzyme_result,
        expression_compatibility=expression_compatibility_result,
        druggability_report=druggability_result,
        conformational_analysis=conformational_result,
        inverse_designs=inverse_designs_result,
        therapeutic_decision=therapeutic_result,
        allosteric_network=allosteric_result,
        epistasis_report=epistasis_result,
        binding_energy_estimates=binding_energy_result or [],
        observables_prediction=observables_result,
        active_learning_recommendation=active_learning_result,
        glycan_analysis=glycan_result,
        admet_report=admet_result,
        directed_evolution_plan=directed_evolution_result,
        protac_feasibility=protac_result,
        selectivity_landscape=selectivity_result,
        fragment_hotspot_map=fragment_hotspot_result,
        aav_design=aav_result,
        crispr_plan=crispr_result,
        genetic_circuit_design=genetic_circuit_result,
        covalent_inhibitor_design=covalent_result,
        engineering_strategy=engineering_result,
        hbond_network=hbond_result,
        disease_association_report=disease_association_result,
        etp_analysis=etp_result,
    )

    # ── AI summary (Gemini) ──
    if ai_summary:
        try:
            from core.gemini_interpreter import generate_report_summary
            record.ai_summary = await generate_report_summary(record)
        except Exception:
            pass

    # ── Download pathway diagrams after record is assembled ──
    if include_pathways and not no_images:
        try:
            from analyzers.diagram_fetcher import fetch_all_diagrams
            await fetch_all_diagrams(
                record,
                save_dir="./diagrams",
                annotate=not no_annotate,
                max_callouts=max_callouts,
                no_ai_polish=no_ai_polish,
            )
        except Exception:
            pass

    return record


async def run_reverse_genetics_query(
    gene: str,
    mutation_type: str,
    organism: str = "Homo sapiens",
    specific_variant: Optional[str] = None,
    step_cb=None,
) -> dict:
    """Run reverse genetics analysis for a gene mutation.

    Returns a dict serialization of ReverseGeneticsReport.
    """
    from analyzers.reverse_genetics_analyzer import run_reverse_genetics_analysis
    report = await run_reverse_genetics_analysis(
        gene=gene,
        mutation_type=mutation_type,
        organism=organism,
        specific_variant=specific_variant,
        step_cb=step_cb,
    )
    return report.model_dump() if hasattr(report, "model_dump") else vars(report)


async def run_rnai_query(
    gene: str,
    rnai_type: str = "siRNA",
    organism: str = "Homo sapiens",
    cell_type: Optional[str] = None,
    step_cb=None,
) -> dict:
    """Run RNAi experiment design analysis for a gene.

    Returns a dict serialization of RNAiReport.
    """
    from analyzers.rnai_analyzer import run_rnai_analysis
    report = await run_rnai_analysis(
        gene=gene,
        rnai_type=rnai_type,
        organism=organism,
        cell_type=cell_type,
        step_cb=step_cb,
    )
    return report.model_dump() if hasattr(report, "model_dump") else vars(report)


async def run_compatibility_query(
    source_gene: str,
    source_organism: str,
    target_organism: str,
    target_cell: Optional[str] = None,
    step_cb=None,
) -> dict:
    """Run cross-species compatibility + metabolic insertion impact analysis.

    Returns a dict with:
        compatibility: CompatibilityReport (serialized)
        insertion_impact: InsertionImpactReport (serialized)
    """
    from analyzers.compatibility_analyzer import run_compatibility_with_impact
    compat, impact = await run_compatibility_with_impact(
        source_gene=source_gene,
        source_organism=source_organism,
        target_organism=target_organism,
        target_cell=target_cell,
        step_cb=step_cb,
    )
    return {
        "compatibility": compat.model_dump() if hasattr(compat, "model_dump") else vars(compat),
        "insertion_impact": impact.model_dump() if hasattr(impact, "model_dump") else vars(impact),
    }


def run_query(
    query: str,
    organism: Optional[str],
    fmt: str,
    output: str,
    cross_species: bool,
    no_cache: bool = False,
    include_pathways: bool = True,
    no_images: bool = False,
    min_pgx_evidence: str = "4",
    interaction_score: float = 0.4,
    ai_summary: bool = False,
    no_annotate: bool = False,
    max_callouts: int = 10,
    no_ai_polish: bool = False,
) -> None:
    """Synchronous entry point called by main.py Click command."""
    import click

    cache = None if no_cache else Cache()

    async def _run():
        accession = await _resolve_accession(query, organism)
        click.echo(f"Resolved accession: {accession}")
        record = await fetch_all(
            accession, organism, cross_species, cache,
            include_pathways=include_pathways,
            no_images=no_images,
            min_pgx_evidence=min_pgx_evidence,
            interaction_score=interaction_score,
            ai_summary=ai_summary,
            no_annotate=no_annotate,
            max_callouts=max_callouts,
            no_ai_polish=no_ai_polish,
        )
        return record

    record = asyncio.run(_run())

    click.echo(f"Generating {fmt} report...")

    if fmt == "json":
        from output.json_writer import write_json
        dest = write_json(record, output)
    elif fmt == "pdf":
        from output.pdf_writer import write_pdf
        dest = write_pdf(record, output)
    else:  # markdown
        from output.markdown_writer import write_markdown
        dest = write_markdown(record, output)

    click.echo(f"Report written to: {dest}")

    if cache:
        cache.close()
