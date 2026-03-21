"""Gemini AI interpreter — NL query resolution and report summarization.

Uses Google Gemini 1.5 Flash (free tier). All functions degrade gracefully
if GEMINI_API_KEY is not set or quota is exceeded — the report still generates
without AI content.

Install: pip install google-generativeai
"""

from __future__ import annotations

import asyncio
import logging
import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

# Load .env from proteinscope/ directory (parent of core/)
load_dotenv(Path(__file__).parent.parent / ".env", override=True)

log = logging.getLogger(__name__)

# Model name (google.genai SDK)
_MODEL_NAME = "models/gemini-2.5-pro"
_client = None
_client_key: str = ""


def _get_client():
    """Lazy-initialize google.genai Client; re-initializes if key changes."""
    global _client, _client_key
    key = os.getenv("GEMINI_API_KEY", "")
    if not key:
        return None
    if _client is not None and key == _client_key:
        return _client
    _client = None
    _client_key = key
    try:
        from google import genai
        _client = genai.Client(api_key=key)
        log.debug("google.genai Client initialized")
        return _client
    except Exception as exc:
        log.warning("google.genai init failed: %s", exc)
    return None


async def _call(prompt: str, timeout: float = 120.0) -> str:
    """Run Gemini 2.5 Pro with a timeout. Returns '' on failure."""
    client = _get_client()
    if client is None:
        return ""
    try:
        loop = asyncio.get_event_loop()
        response = await asyncio.wait_for(
            loop.run_in_executor(
                None,
                lambda: client.models.generate_content(
                    model=_MODEL_NAME, contents=prompt
                ),
            ),
            timeout=timeout,
        )
        return (response.text or "").strip()
    except asyncio.TimeoutError:
        log.warning("Gemini timed out after %.0fs", timeout)
        return ""
    except Exception as exc:
        log.warning("Gemini call failed: %s", exc)
        return ""


# Utility calls (translate/classify/parse) use same model, shorter timeout
async def _call_fast(prompt: str) -> str:
    return await _call(prompt, timeout=45.0)


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

async def classify_query_intent(query: str, *, session_context: str = "") -> dict:
    """Top-level intent classification — runs BEFORE any protein resolution.

    Determines whether the user's query is something ProteinScope can process
    or whether it is off-topic / lacks the biological context needed.

    Returns a dict with:
        intent (str): one of
            "protein_query"   — biological/protein question ProteinScope can handle
            "off_topic"       — no protein/biological context; show usage guide
        reason (str): short explanation (for logging / debug)

    Degrades gracefully: returns {"intent": "protein_query", "reason": "fallback"}
    if Gemini is unavailable (so the existing pipeline handles it).
    """
    import json as _json
    prompt = (
        "You are the intent classifier for ProteinScope, a protein analysis platform.\n\n"
        "ProteinScope CAN handle:\n"
        "  - Queries about specific proteins or genes (e.g. 'BRCA1', 'insulin receptor')\n"
        "  - Biological process/pathway queries (e.g. 'DNA replication proteins in penguin')\n"
        "  - Cross-species compatibility questions (e.g. 'Can penguin POLD1 work in human cells?')\n"
        "  - Protein report or analysis requests WITH a protein/gene/pathway mentioned\n"
        "  - Any query that contains a protein name, gene symbol, enzyme, pathway, or organism\n"
        "  - Queries in ANY language (Korean, Japanese, Chinese, Hebrew, Arabic, etc.)\n\n"
        "ProteinScope CANNOT handle:\n"
        "  - Requests that reference 'this' or 'it' without naming a protein or pathway\n"
        "  - Vague follow-up commands with no biological entity (e.g. 'Make a report on this', "
        "'Why is this impossible?', 'Explain that')\n"
        "  - General knowledge questions unrelated to proteins (weather, jokes, math, etc.)\n\n"
        f"User query: \"{query}\"\n\n"
        "Classify this query. Return ONLY a raw JSON object — no markdown, no explanation:\n"
        "  {\"intent\": \"protein_query\", \"reason\": \"<1-sentence reason>\"}\n"
        "  OR\n"
        "  {\"intent\": \"off_topic\", \"reason\": \"<1-sentence reason>\"}\n"
    )
    if session_context:
        prompt = session_context + "\n\n" + prompt
    result = await _call_fast(prompt)
    if not result:
        return {"intent": "protein_query", "reason": "fallback — Gemini unavailable"}
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        intent = str(data.get("intent", "protein_query")).strip()
        if intent not in ("protein_query", "off_topic"):
            intent = "protein_query"
        return {"intent": intent, "reason": str(data.get("reason", ""))}
    except Exception:
        return {"intent": "protein_query", "reason": "fallback — parse error"}


async def translate_to_english(user_query: str) -> str:
    """Translate a query to English if it is written in another language.

    Detects non-English text by counting non-ASCII characters.
    If the query is already English (< 20% non-ASCII), returns it unchanged.
    Falls back to original query if translation fails.
    """
    non_ascii = sum(1 for c in user_query if ord(c) > 127)
    if non_ascii < len(user_query) * 0.15:
        return user_query  # already English / mostly ASCII

    prompt = (
        "Translate the following query into English. "
        "The query may be in any language (Korean, Japanese, Chinese, Hebrew, Arabic, "
        "Russian, Spanish, French, German, etc.). "
        "Return ONLY the English translation — no explanation, no quotes:\n\n"
        f"{user_query}"
    )
    result = await _call_fast(prompt)
    return result.strip() if result else user_query


async def parse_nl_protein_query(user_query: str, *, session_context: str = "") -> dict:
    """Parse a natural language protein query to understand intent and context.

    Returns a dict with:
        is_multi (bool): True if the query is about a process/pathway/multiple proteins.
        gene (str|None): Gene symbol if the query clearly refers to one specific protein.
        organism (str|None): Organism mentioned in the query (e.g. "penguin", "mouse").

    Examples:
        "enzyme that breaks down warfarin"          -> {"is_multi": False, "gene": "CYP2C9", "organism": None}
        "all proteins in DNA replication"           -> {"is_multi": True,  "gene": None,    "organism": None}
        "proteins in DNA replication in penguin"    -> {"is_multi": True,  "gene": None,    "organism": "penguin"}
        "BRCA1 in mouse"                            -> {"is_multi": False, "gene": "BRCA1", "organism": "mouse"}
    """
    import json as _json
    prompt = (
        "You are a molecular biology expert parsing a protein database search query.\n"
        "The query may be written in ANY language — English, Korean, Japanese, Chinese, "
        "Hebrew, Arabic, Russian, Spanish, French, German, or any other language. "
        "Understand the semantic meaning regardless of language or script.\n\n"
        f"User query: \"{user_query}\"\n\n"
        "Analyze this query and return a JSON object with exactly these fields:\n"
        "  is_multi: true if the query mentions a biological process, pathway, or multiple proteins "
        "(e.g. 'DNA replication', 'glycolysis', 'all proteins in apoptosis', "
        "'사과향을 만드는 단백질들', 'アポトーシスのタンパク質', 'חלבוני שכפול ה-DNA', 'بروتينات تكرار الحمض النووي'); "
        "false if it refers to a single specific protein or gene.\n"
        "  gene: the gene symbol (e.g. CYP2C9, EGFR) if is_multi is false and you are confident; "
        "null otherwise.\n"
        "  organism: the organism name in English if explicitly mentioned in ANY language "
        "(e.g. 'penguin' from 팽귄/ペンギン/企鹅/פינגווין, 'mouse', 'apple tree'); null if not specified.\n\n"
        "Return ONLY the raw JSON object — no markdown, no explanation:\n"
        '{"is_multi": true, "gene": null, "organism": "penguin"}'
    )
    if session_context:
        prompt = session_context + "\n\n" + prompt
    result = await _call_fast(prompt)
    if not result:
        return {"is_multi": False, "gene": None, "organism": None}
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        return {
            "is_multi": bool(data.get("is_multi", False)),
            "gene": str(data["gene"]).upper().strip() if data.get("gene") else None,
            "organism": str(data["organism"]).strip() if data.get("organism") else None,
        }
    except Exception:
        return {"is_multi": False, "gene": None, "organism": None}


async def resolve_nl_query(user_query: str) -> Optional[str]:
    """Convert a natural language protein description to a gene symbol or accession.

    Returns a gene symbol / UniProt accession string, or None if not confident.

    Examples:
        "the enzyme that breaks down warfarin"  -> "CYP2C9"
        "cystic fibrosis channel protein"       -> "CFTR"
        "HER2 cancer target"                    -> "ERBB2"
        "all proteins in DNA replication"       -> "PCNA"  (most representative)
    """
    prompt = (
        "You are a molecular biology expert helping a protein database search tool.\n"
        "The tool searches ONE protein at a time from UniProt.\n"
        "The query may be in ANY language (English, Korean, Japanese, Chinese, etc.). "
        "Translate and interpret it, then convert to a SINGLE gene symbol.\n\n"
        "Rules:\n"
        "- Return ONLY the gene symbol (e.g. PCNA, EGFR, TP53) — no explanation, no spaces, no punctuation.\n"
        "- If the query describes a biological process or pathway (e.g. 'DNA replication', "
        "'glycolysis', 'apoptosis', '사과향', 'apple scent biosynthesis'), return the MOST "
        "REPRESENTATIVE or MOST STUDIED single gene symbol for that process.\n"
        "- For scent/aroma biosynthesis queries, return the key biosynthetic enzyme "
        "(e.g. apple scent → MdAAT1, rose scent → RcRAS1).\n"
        "- If the query is completely unrelated to biology or proteins, reply with exactly: UNKNOWN\n\n"
        f"Description: {user_query}"
    )
    result = await _call_fast(prompt)
    if not result or result.upper() == "UNKNOWN":
        return None
    # Strip any accidental whitespace or punctuation; take first token only
    token = result.split()[0].strip(".,;:'\"()[]").upper()
    return token if token and token != "UNKNOWN" else None


async def suggest_proteins_for_query(
    user_query: str,
    organism: Optional[str] = None,
) -> list[dict]:
    """Return a list of {gene, description} dicts for proteins related to the query.

    Used when a query cannot be resolved to a single UniProt entry to give the
    user clickable protein options in the web UI.

    Args:
        user_query: The original user query string.
        organism:   Organism context extracted from the query (e.g. "penguin").
                    If provided, suggestions will note whether the organism has
                    limited UniProt coverage and show the best available orthologs.

    Returns a list like:
        [{"gene": "PCNA", "description": "Sliding clamp / processivity factor"}, ...]
    Returns [] on failure.
    """
    import json as _json
    org_context = (
        f"The user is interested in the organism: {organism}. "
        "List gene symbols using their standard names (e.g. PCNA, POLD1) — do NOT switch to human. "
        "The search system will look for each gene in the specified organism first and fall back "
        "automatically if that organism's entry is not available in UniProt.\n"
    ) if organism else ""
    prompt = (
        "You are a molecular biology expert helping a user search a protein database.\n"
        "The user's query matches multiple proteins or a biological process.\n"
        "The query may be in ANY language — understand it and respond with standard gene symbols.\n\n"
        f"User query: \"{user_query}\"\n"
        f"{org_context}\n"
        "Return a JSON array of 5 key proteins with their standard gene symbols and one-line "
        "functional descriptions in English. Format EXACTLY like this — no markdown fences, no explanation, "
        "ONLY the raw JSON array:\n"
        '[{"gene": "PCNA", "description": "Sliding clamp / DNA polymerase processivity factor"}, ...]\n\n'
        "Use real, well-known gene symbols that exist in UniProt."
    )
    result = await _call_fast(prompt)
    if not result:
        return []
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        if isinstance(data, list):
            return [
                {"gene": str(d.get("gene", "")).upper().strip(),
                 "description": str(d.get("description", "")).strip()}
                for d in data
                if d.get("gene")
            ]
    except Exception:
        pass
    return []


async def generate_report_summary(record) -> str:
    """Generate a plain-English summary of a ProteinRecord (max ~300 words).

    Covers: biological function, clinical significance, key drugs,
    pathway roles, and notable research highlights.
    Returns empty string if key is missing or generation fails.
    """
    diseases = ", ".join(
        list({cv.disease for cv in record.clinical_variants})[:5]
    ) or "none recorded"

    drug_names = ", ".join(
        [di.drug_name for di in record.drug_interactions[:8]]
    ) or "none recorded"

    pathway_names = ", ".join(
        [mp.pathway_name for mp in (record.metabolic_pathways + record.signaling_pathways)[:6]]
    ) or "none recorded"

    pgx_summary = (
        f"{len(record.pgx_variants)} PGx variant(s) recorded"
        if record.pgx_variants
        else "none recorded"
    )

    prompt = (
        "You are a biomedical writer summarizing a protein research report.\n"
        "Write a concise plain-English summary (max 300 words) covering:\n"
        "1. What this protein does biologically\n"
        "2. Its clinical significance and associated diseases\n"
        "3. Key drug interactions or pharmacogenomic variants (if any)\n"
        "4. Its role in major signaling or metabolic pathways (if any)\n"
        "5. Any notable research highlights\n\n"
        f"Protein: {record.protein_name} ({record.gene_name}), {record.organism}\n"
        f"Function: {record.function_description[:500]}\n"
        f"Diseases: {diseases}\n"
        f"Key drugs: {drug_names}\n"
        f"Pathways: {pathway_names}\n"
        f"PGx variants: {pgx_summary}\n\n"
        "Write in a factual, academic tone. Do not invent data not listed above."
    )
    return await _call(prompt)


async def parse_compatibility_query(user_query: str) -> dict:
    """Detect if a query asks about cross-species protein substitution.

    Returns dict with:
        is_compatibility: bool
        source_gene: str | None       (gene to transplant)
        source_organism: str | None   (where it comes from)
        target_organism: str | None   (where it goes)
        target_cell: str | None       (specific cell line, e.g. "HeLa")

    Examples:
        "Can I use penguin POLD1 in HeLa cells?"
            -> {is_compatibility: True, source_gene: "POLD1", source_organism: "penguin",
               target_organism: "Homo sapiens", target_cell: "HeLa"}
        "Would mouse EGFR work in human cells?"
            -> {is_compatibility: True, source_gene: "EGFR", source_organism: "mouse",
               target_organism: "Homo sapiens", target_cell: None}
        "What is EGFR?"
            -> {is_compatibility: False, ...}
    """
    import json as _json
    prompt = (
        "You are a molecular biology expert parsing a protein substitution or heterologous expression query.\n"
        "The query may be in ANY language (English, Korean, Japanese, Chinese, etc.). "
        "Understand the semantic meaning regardless of language.\n\n"
        f"User query: \"{user_query}\"\n\n"
        "Determine if this query asks whether a protein from one species/organism can be expressed "
        "or functionally work in cells of another species/organism. This includes:\n"
        "- Cross-species substitution: 'Can penguin POLD1 replace human POLD1 in HeLa cells?'\n"
        "- Heterologous expression: 'Can apple scent enzyme MdAAT1 be expressed in human cells?'\n"
        "- Plant→animal: 'Can I express a plant enzyme in human/animal cells?'\n\n"
        "Return ONLY a raw JSON object with these fields:\n"
        "  is_compatibility: true if this is a cross-species expression or substitution question\n"
        "  source_gene: the gene symbol in English (e.g. POLD1, MdAAT1, EGFR). "
        "Resolve common or non-English descriptions to the best known gene symbol: "
        "'사과향 단백질' or 'apple scent protein' → 'MdAAT1'; "
        "'DNA polymerase' → 'POLD1'; 'rose scent' → 'RcRAS1'. "
        "Null only if truly unresolvable.\n"
        "  source_organism: organism the protein comes FROM, in English (e.g. 'Malus domestica', 'penguin', 'apple tree'); else null\n"
        "  target_organism: organism/system to express the protein IN, in English (e.g. 'Homo sapiens', 'human'); else null\n"
        "  target_cell: specific cell line if mentioned (e.g. 'HeLa', 'HEK293'); else null\n\n"
        "Examples:\n"
        '{"is_compatibility": true, "source_gene": "POLD1", "source_organism": "penguin", '
        '"target_organism": "Homo sapiens", "target_cell": "HeLa"}\n'
        '{"is_compatibility": true, "source_gene": "MdAAT1", "source_organism": "Malus domestica", '
        '"target_organism": "Homo sapiens", "target_cell": null}\n\n'
        "Return ONLY the JSON, no markdown:"
    )
    result = await _call_fast(prompt)
    if not result:
        return {"is_compatibility": False, "source_gene": None, "source_organism": None,
                "target_organism": None, "target_cell": None}
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        return {
            "is_compatibility": bool(data.get("is_compatibility", False)),
            "source_gene": str(data["source_gene"]).upper().strip() if data.get("source_gene") else None,
            "source_organism": str(data["source_organism"]).strip() if data.get("source_organism") else None,
            "target_organism": str(data["target_organism"]).strip() if data.get("target_organism") else None,
            "target_cell": str(data["target_cell"]).strip() if data.get("target_cell") else None,
        }
    except Exception:
        return {"is_compatibility": False, "source_gene": None, "source_organism": None,
                "target_organism": None, "target_cell": None}


async def parse_reverse_genetics_query(user_query: str) -> dict:
    """Detect if a query asks about the effect of a gene mutation (reverse genetics).

    Returns:
        is_reverse_genetics: bool
        gene: str | None           gene symbol
        mutation_type: str         loss_of_function | gain_of_function | missense |
                                   dominant_negative | haploinsufficiency | frameshift |
                                   nonsense | deletion | amplification | general
        specific_variant: str|None HGVS notation if mentioned (e.g. "p.R175H")
        organism: str|None
    """
    import json as _json
    prompt = (
        "You are a molecular genetics expert parsing a reverse genetics query.\n"
        "The query may be in ANY language. Detect the semantic meaning.\n\n"
        f"User query: \"{user_query}\"\n\n"
        "Determine if this query asks: 'What happens when gene X is mutated/knocked out/deleted/amplified?'\n"
        "This includes:\n"
        "  - Knockout/knockdown: 'What if TP53 is knocked out?', 'TP53 결손이 생기면?'\n"
        "  - Specific mutation: 'BRCA1 R175H mutation analysis', 'What does EGFR T790M do?'\n"
        "  - LoF: 'What happens if PTEN is inactivated?'\n"
        "  - GoF: 'KRAS G12D gain of function effects'\n"
        "  - General: 'What pathways are affected by ATM mutation?'\n\n"
        "Return ONLY raw JSON:\n"
        "  is_reverse_genetics: true/false\n"
        "  gene: gene symbol in English (e.g. BRCA1, TP53) or null\n"
        "  mutation_type: one of loss_of_function | gain_of_function | missense | "
        "dominant_negative | haploinsufficiency | frameshift | nonsense | deletion | amplification | general\n"
        "  specific_variant: HGVS string if mentioned (e.g. 'p.R175H') else null\n"
        "  organism: organism in English else null\n\n"
        "Examples:\n"
        '{"is_reverse_genetics":true,"gene":"TP53","mutation_type":"loss_of_function","specific_variant":null,"organism":"Homo sapiens"}\n'
        '{"is_reverse_genetics":true,"gene":"BRCA1","mutation_type":"missense","specific_variant":"p.R175H","organism":null}\n'
        '{"is_reverse_genetics":false,"gene":null,"mutation_type":"general","specific_variant":null,"organism":null}\n\n'
        "Return ONLY the JSON, no markdown:"
    )
    result = await _call_fast(prompt)
    _default = {"is_reverse_genetics": False, "gene": None, "mutation_type": "general",
                "specific_variant": None, "organism": None}
    if not result:
        return _default
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        valid_mt = {
            "loss_of_function", "gain_of_function", "missense", "dominant_negative",
            "haploinsufficiency", "frameshift", "nonsense", "deletion", "amplification", "general",
        }
        mt = str(data.get("mutation_type", "general")).strip().lower()
        return {
            "is_reverse_genetics": bool(data.get("is_reverse_genetics", False)),
            "gene": str(data["gene"]).upper().strip() if data.get("gene") else None,
            "mutation_type": mt if mt in valid_mt else "general",
            "specific_variant": str(data["specific_variant"]).strip() if data.get("specific_variant") else None,
            "organism": str(data["organism"]).strip() if data.get("organism") else None,
        }
    except Exception:
        return _default


async def parse_rnai_query(user_query: str) -> dict:
    """Detect if a query requests an RNAi knockdown experiment design.

    Returns:
        is_rnai: bool
        gene: str | None           target gene symbol
        rnai_type: str             siRNA | shRNA | miRNA | general
        cell_type: str | None      target cell type/line
        organism: str | None
    """
    import json as _json
    prompt = (
        "You are a molecular biology expert parsing an RNAi experiment design query.\n"
        "The query may be in ANY language.\n\n"
        f"User query: \"{user_query}\"\n\n"
        "Determine if this query asks to design, analyze, or interpret an RNAi (RNA interference) "
        "experiment — siRNA, shRNA, or miRNA knockdown of a specific gene.\n"
        "This includes:\n"
        "  - 'Design siRNA for BRCA1 knockdown in HeLa cells'\n"
        "  - 'EGFR RNAi experiment in lung cancer cells'\n"
        "  - 'TP53 shRNA knockdown — what would happen?'\n"
        "  - 'BRCA1 RNAi 설계해줘'\n\n"
        "Return ONLY raw JSON:\n"
        "  is_rnai: true/false\n"
        "  gene: gene symbol in English or null\n"
        "  rnai_type: siRNA | shRNA | miRNA | general\n"
        "  cell_type: cell line or tissue (e.g. 'HeLa', 'lung cancer', 'HEK293') or null\n"
        "  organism: organism in English or null\n\n"
        "Return ONLY the JSON, no markdown:"
    )
    result = await _call_fast(prompt)
    _default = {"is_rnai": False, "gene": None, "rnai_type": "general",
                "cell_type": None, "organism": None}
    if not result:
        return _default
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        valid_type = {"siRNA", "shRNA", "miRNA", "general"}
        rt = str(data.get("rnai_type", "general")).strip()
        return {
            "is_rnai": bool(data.get("is_rnai", False)),
            "gene": str(data["gene"]).upper().strip() if data.get("gene") else None,
            "rnai_type": rt if rt in valid_type else "general",
            "cell_type": str(data["cell_type"]).strip() if data.get("cell_type") else None,
            "organism": str(data["organism"]).strip() if data.get("organism") else None,
        }
    except Exception:
        return _default


async def extract_organism_from_query(original_query: str, english_query: str) -> Optional[str]:
    """Extract organism name from a query, trying both original and English forms.

    Returns the organism name in English (e.g. "penguin", "apple tree") or None.
    This is a focused single-purpose call, more reliable than parsing organism
    as a side-effect of parse_nl_protein_query.
    """
    prompt = (
        "You are a biologist. The user asked a protein database question.\n"
        f"Original query: \"{original_query}\"\n"
        f"English translation: \"{english_query}\"\n\n"
        "Is a specific organism (species or taxonomic group) mentioned in either version of the query?\n"
        "If yes, return ONLY the organism name in English (e.g. 'penguin', 'apple tree', 'mouse', 'E. coli').\n"
        "If no organism is mentioned, return: null\n\n"
        "Return ONLY the organism name or null — no explanation, no quotes:\n"
        "penguin"
    )
    result = await _call_fast(prompt)
    if not result:
        return None
    cleaned = result.strip().strip('"').strip("'").lower()
    if cleaned in ("null", "none", "", "no", "not specified", "not mentioned"):
        return None
    return cleaned.capitalize() if cleaned else None


async def detect_taxon_ambiguity(organism: str) -> dict:
    """Detect if an organism name refers to a broad taxonomic group.

    Returns dict with:
        is_broad: bool — True if genus/family/order/class level
        taxon_rank: str — "species"|"genus"|"family"|"order"|"class"
        species: list[dict] — up to 6 notable species if is_broad is True
            each: {common_name, scientific_name, note}
    """
    import json as _json
    prompt = (
        f"You are a taxonomy expert. The user mentioned organism: '{organism}'\n\n"
        "Determine if this refers to a BROAD taxonomic group (genus, family, order, class) "
        "rather than a single specific species.\n\n"
        "Broad examples: 'penguin' (order Sphenisciformes), 'snake' (suborder Serpentes), "
        "'fish' (superclass), 'bacteria', 'pine tree' (genus Pinus), 'oak'.\n"
        "Specific examples: 'mouse' (Mus musculus), 'rat' (Rattus norvegicus), "
        "'human' (Homo sapiens), 'E. coli', 'yeast' (Saccharomyces cerevisiae), "
        "'zebrafish' (Danio rerio) — these are conventionally treated as single model organisms.\n\n"
        "If IS broad: list the 5-6 most studied / notable species in that group.\n"
        "If NOT broad: return is_broad: false with no species list.\n\n"
        "Return ONLY raw JSON — no markdown, no explanation:\n"
        '{"is_broad": true, "taxon_rank": "order", "species": ['
        '{"common_name": "Emperor Penguin", "scientific_name": "Aptenodytes forsteri", "note": "Largest penguin; most studied for cold-adaptation proteins"},'
        '{"common_name": "Adélie Penguin", "scientific_name": "Pygoscelis adeliae", "note": "Model for Antarctic ecology; sequenced genome"}'
        ']}\n'
        'OR\n{"is_broad": false, "taxon_rank": "species", "species": []}'
    )
    result = await _call_fast(prompt)
    if not result:
        return {"is_broad": False, "taxon_rank": "species", "species": []}
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        return {
            "is_broad": bool(data.get("is_broad", False)),
            "taxon_rank": str(data.get("taxon_rank", "species")),
            "species": [
                {
                    "common_name": str(s.get("common_name", "")),
                    "scientific_name": str(s.get("scientific_name", "")),
                    "note": str(s.get("note", "")),
                }
                for s in data.get("species", [])
                if s.get("scientific_name")
            ],
        }
    except Exception:
        return {"is_broad": False, "taxon_rank": "species", "species": []}


async def parse_pathway_analysis_plan(user_query: str) -> dict:
    """Detect if a query asks to engineer/implement a complex multi-protein biological pathway.

    TRIGGERS when ALL three are true:
    1. User wants to IMPLEMENT/ENABLE a complete biological pathway (not just analyze one protein)
    2. The pathway requires multiple interdependent proteins working together
    3. The proteins are foreign to the target organism/cell

    Does NOT trigger on:
    - Single protein compatibility ("Can RBCL work in human cells?")
    - Pathway information queries ("what proteins are in glycolysis?")
    - RNAi / knockdown / mutation analysis

    Returns dict with:
        is_pathway_analysis: bool
        pathway_name:        str
        overview:            str  (1-2 sentence engineering challenge summary)
        target_organism:     str
        target_cell:         str | null
        source_organism:     str
        analysis_steps:      list[dict]  (4-8 steps)
    """
    import json as _json

    _default = {
        "is_pathway_analysis": False, "pathway_name": "", "overview": "",
        "target_organism": "Homo sapiens", "target_cell": None, "source_organism": "",
        "analysis_steps": [],
    }

    prompt = (
        "You are a systems biology expert. Determine if this query asks to ENGINEER or IMPLEMENT "
        "a complete multi-protein biological pathway in a foreign cellular context.\n\n"
        "TRIGGER (ALL must be true):\n"
        "1. User wants to implement/enable/engineer a complete biological pathway\n"
        "2. That pathway requires multiple proteins working together\n"
        "3. The proteins come from a different organism or are otherwise foreign to the target cell\n\n"
        "DOES NOT TRIGGER:\n"
        "- Single protein compatibility ('Can RBCL work in human cells?')\n"
        "- Pathway information queries ('what are the glycolysis enzymes?')\n"
        "- RNAi, knockdown, or mutation analysis\n\n"
        "TRIGGERS:\n"
        "- 'Enable photosynthesis in human epithelial cells'\n"
        "- 'Engineer nitrogen fixation in yeast'\n"
        "- 'Express the complete mevalonate pathway in E. coli'\n"
        "- '사람 상피세포에서 광합성을 할 수 있도록 하게 광합성에 관련된것을 사람 상피세포에서 발현시키고 싶어'\n\n"
        f"Query: \"{user_query}\"\n\n"
        "If is_pathway_analysis is TRUE, generate 4-8 analysis steps:\n"
        "- Step 1: ALWAYS type='info' — write 200-word biological overview listing ALL proteins needed, "
        "organized by functional module (e.g. Photosystem II proteins, PSI proteins, Calvin cycle, etc.)\n"
        "- Middle steps: type='compatibility' — pick 3-5 REPRESENTATIVE proteins/complexes to analyze "
        "(most important ones; don't repeat ones already checked)\n"
        "- Last step: ALWAYS type='summary' — overall feasibility synthesis\n\n"
        "Return ONLY raw JSON:\n"
        "{\n"
        '  "is_pathway_analysis": true,\n'
        '  "pathway_name": "Photosynthesis Engineering in Human Epithelial Cells",\n'
        '  "overview": "1-2 sentence description of the engineering challenge",\n'
        '  "target_organism": "Homo sapiens",\n'
        '  "target_cell": "epithelial cell",\n'
        '  "source_organism": "Spinacia oleracea / Arabidopsis thaliana",\n'
        '  "analysis_steps": [\n'
        '    {"id": 1, "title": "Pathway Component Map", "description": "Overview of all required components", '
        '"type": "info", "content": "200-word biological overview..."},\n'
        '    {"id": 2, "title": "Photosystem II Core (PsbA/D1)", "description": "PSII reaction center", '
        '"type": "compatibility", "params": {"source_gene": "PSBA", "source_organism": "Spinacia oleracea", '
        '"target_organism": "Homo sapiens", "target_cell": "epithelial cell"}},\n'
        '    {"id": 7, "title": "Overall Feasibility Assessment", '
        '"description": "Synthesis of all analyses", "type": "summary"}\n'
        "  ]\n"
        "}"
    )

    result = await _call_fast(prompt)
    if not result:
        return _default
    try:
        cleaned = result.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        if not data.get("is_pathway_analysis"):
            return _default

        steps = []
        for s in data.get("analysis_steps", []):
            step_type = s.get("type", "info")
            if step_type not in ("info", "compatibility", "summary"):
                step_type = "info"
            step = {
                "id": int(s.get("id", len(steps) + 1)),
                "title": str(s.get("title", "")),
                "description": str(s.get("description", "")),
                "type": step_type,
            }
            if step_type == "info":
                step["content"] = str(s.get("content", ""))
            elif step_type == "compatibility":
                step["params"] = {
                    "source_gene": str(s.get("params", {}).get("source_gene", "")),
                    "source_organism": str(s.get("params", {}).get("source_organism", "")),
                    "target_organism": str(s.get("params", {}).get("target_organism", "")),
                    "target_cell": s.get("params", {}).get("target_cell") or None,
                }
            steps.append(step)

        return {
            "is_pathway_analysis": True,
            "pathway_name": str(data.get("pathway_name", "")),
            "overview": str(data.get("overview", "")),
            "target_organism": str(data.get("target_organism", "Homo sapiens")),
            "target_cell": data.get("target_cell") or None,
            "source_organism": str(data.get("source_organism", "")),
            "analysis_steps": steps,
            "references": _PATHWAY_REFERENCES,
        }
    except Exception:
        return _default


_PATHWAY_REFERENCES: list[dict] = [
    {
        "pubmed_id": "26967285",
        "doi": "10.1016/j.cell.2016.02.004",
        "title": "Engineering cellular metabolism",
        "authors": ["Nielsen J", "Keasling JD"],
        "journal": "Cell",
        "year": 2016,
        "apa_citation": "Nielsen, J., & Keasling, J. D. (2016). Engineering cellular metabolism. Cell, 164(6), 1185–1197. https://doi.org/10.1016/j.cell.2016.02.004",
    },
    {
        "pubmed_id": "16612385",
        "doi": "10.1038/nature04640",
        "title": "Production of the antimalarial drug precursor artemisinic acid in engineered yeast",
        "authors": ["Ro DK", "Paradise EM", "Ouellet M", "Fisher KJ", "Newman KL", "Ndungu JM", "Keasling JD"],
        "journal": "Nature",
        "year": 2006,
        "apa_citation": "Ro, D. K., Paradise, E. M., Ouellet, M., Fisher, K. J., Newman, K. L., Ndungu, J. M., & Keasling, J. D. (2006). Production of the antimalarial drug precursor artemisinic acid in engineered yeast. Nature, 440(7086), 940–943. https://doi.org/10.1038/nature04640",
    },
    {
        "pubmed_id": "22745249",
        "doi": "10.1126/science.1225829",
        "title": "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity",
        "authors": ["Jinek M", "Chylinski K", "Fonfara I", "Hauer M", "Doudna JA", "Charpentier E"],
        "journal": "Science",
        "year": 2012,
        "apa_citation": "Jinek, M., Chylinski, K., Fonfara, I., Hauer M., Doudna, J. A., & Charpentier, E. (2012). A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. Science, 337(6096), 816–821. https://doi.org/10.1126/science.1225829",
    },
    {
        "pubmed_id": "22342777",
        "doi": "10.1016/j.ymben.2011.12.004",
        "title": "Parts plus pipes: Synthetic biology approaches to metabolic engineering",
        "authors": ["Boyle PM", "Silver PA"],
        "journal": "Metabolic Engineering",
        "year": 2012,
        "apa_citation": "Boyle, P. M., & Silver, P. A. (2012). Parts plus pipes: Synthetic biology approaches to metabolic engineering. Metabolic Engineering, 14(3), 223–232. https://doi.org/10.1016/j.ymben.2011.12.004",
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


async def generate_pathway_summary(plan: dict, completed_steps: list[dict]) -> str:
    """Generate a systems-biology synthesis of all pathway engineering step results.

    Args:
        plan:            The original pathway plan dict.
        completed_steps: List of {id, title, result} dicts from completed analyses.

    Returns:
        ~300-word scientific feasibility assessment string.
    """
    context_lines = [
        f"Pathway Engineering Goal: {plan.get('pathway_name', '')}",
        f"Overview: {plan.get('overview', '')}",
        f"Source organism: {plan.get('source_organism', '')}",
        f"Target organism/cell: {plan.get('target_organism', '')} — {plan.get('target_cell') or 'unspecified cell type'}",
        "",
        "Step Results:",
    ]
    for s in completed_steps:
        context_lines.append(f"\nStep {s['id']}: {s['title']}")
        result = s.get("result", {})
        compat = result.get("compatibility", result)  # handle both flat and nested
        if compat:
            context_lines.append(f"  Overall: {compat.get('overall_assessment', compat.get('gemini_analysis', ''))[:300]}")
            context_lines.append(f"  Compatible: {compat.get('is_compatible', 'unknown')}")
            caveats = compat.get("caveats", [])
            if caveats:
                context_lines.append(f"  Key issues: {'; '.join(str(c) for c in caveats[:2])}")
        impact = result.get("insertion_impact", {})
        if impact and impact.get("overall_metabolic_risk"):
            context_lines.append(f"  Metabolic risk: {impact['overall_metabolic_risk']}")

    context = "\n".join(context_lines)

    prompt = (
        "You are a systems biology and synthetic biology expert. Synthesize the following "
        "pathway engineering analysis into a comprehensive feasibility assessment.\n\n"
        f"{context}\n\n"
        "Provide a rigorous ~300-word assessment covering:\n"
        "1. Overall feasibility verdict (one of: Feasible | Partially Feasible | "
        "Not Feasible in Current State | Theoretically Conceivable but Technically Extreme)\n"
        "2. Principal scientific obstacles (top 3-4 bottlenecks from the data)\n"
        "3. Minimum prerequisites for any progress\n"
        "4. Current state of the art (any published attempts)\n"
        "5. Recommended research priority order\n\n"
        "Write in rigorous scientific style. Base conclusions on the step results provided."
    )
    return await _call(prompt)


async def generate_section_insight(section: str, data_summary: str) -> str:
    """Generate a 2-3 sentence clinical/research insight for a report section.

    Args:
        section:      Section name, e.g. "Drug Interactions", "Disease Cascade".
        data_summary: Brief text description of the data in that section.

    Returns:
        2-3 sentence insight string, or empty string on failure.
    """
    prompt = (
        f"You are a biomedical expert writing a brief research insight.\n"
        f"Section: {section}\n"
        f"Data: {data_summary}\n\n"
        "Write 2-3 concise sentences interpreting the clinical or research "
        "significance of this data. Be factual and academic. "
        "Do not invent information not present in the data above."
    )
    return await _call(prompt)


async def synthesize_ptm_logic(sites_json: str, gene: str) -> dict:
    """Ask Gemini to identify PTM regulatory relationships for a gene.

    Args:
        sites_json: JSON string of PTM site dicts (from fetch_ptm_sites).
        gene:       Gene symbol (e.g. "AKT1").

    Returns a dict with keys:
        state_transitions (list[dict]):  [{from_site, to_site, relation, confidence, evidence}]
        regulatory_summary (str):        Narrative summary of the PTM regulatory logic.
        key_regulatory_nodes (list[str]): Site IDs that are master regulators.
        therapeutic_targets (list[str]): Actionable therapeutic implications.

    Returns empty dict on failure (Gemini unavailable, parse error, etc.).
    """
    import json as _json

    prompt = (
        f"You are a protein biochemistry expert. Given these PTM sites for {gene}:\n"
        f"{sites_json}\n\n"
        "Identify regulatory relationships between sites. Return JSON:\n"
        "{\n"
        '  "state_transitions": [\n'
        '    {"from_site": "T308", "to_site": "S473", "relation": "requires_prior", '
        '"confidence": 0.9, '
        '"evidence": "AKT activation requires PDK1 phosphorylation first"}\n'
        "  ],\n"
        '  "regulatory_summary": "...",\n'
        '  "key_regulatory_nodes": ["T308", "S473"],\n'
        '  "therapeutic_targets": ["T308 — PDK1 inhibition blocks AKT activation cascade"]\n'
        "}\n\n"
        "Return ONLY raw JSON — no markdown fences, no explanation."
    )

    raw = await _call(prompt)
    if not raw:
        return {}

    try:
        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        return {
            "state_transitions": data.get("state_transitions") or [],
            "regulatory_summary": str(data.get("regulatory_summary", "")).strip(),
            "key_regulatory_nodes": [
                str(n).strip() for n in (data.get("key_regulatory_nodes") or []) if n
            ],
            "therapeutic_targets": [
                str(t).strip() for t in (data.get("therapeutic_targets") or []) if t
            ],
        }
    except Exception as exc:
        log.debug("synthesize_ptm_logic parse failed for %s: %s", gene, exc)
        return {}


async def synthesize_therapeutic_strategy(modality_scores: dict, gene: str) -> dict:
    """Gemini synthesis for therapeutic modality comparison (P10)."""
    import json
    try:
        prompt = (
            f"Gene: {gene}\n"
            f"Modality scores: {json.dumps(modality_scores, indent=2)}\n\n"
            "You are a clinical-stage drug developer. Based on this target biology data:\n"
            "1. Which therapeutic modality has the highest probability of clinical success?\n"
            "2. What are the key evidence gaps that could de-risk development?\n"
            "3. Propose the next 3 experiments to validate the lead modality.\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{\"recommended_modality\": \"...\", \"rationale\": \"...\", '
            '\"evidence_gaps\": [\"...\"], \"next_experiments\": [\"...\"]}'
        )
        raw = await _call(prompt)
        if not raw:
            return {}
        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        return json.loads(cleaned)
    except Exception:
        return {}

async def parse_antigen_discovery_query(text: str) -> dict:
    """Extract tumor_type and modality from natural language antigen discovery query (P11)."""
    import json
    try:
        prompt = (
            f"Query: {text}\n\n"
            "Extract the cancer/tumor type and therapeutic modality from this query.\n"
            "Common modalities: CAR-T, CAR-NK, TCE, ADC, antibody.\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{\"tumor_type\": \"melanoma\", \"modality\": \"CAR-T\"}'
        )
        raw = await _call(prompt)
        if not raw:
            return {"tumor_type": "unknown", "modality": "CAR-T"}
        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        return json.loads(cleaned)
    except Exception:
        return {"tumor_type": "unknown", "modality": "CAR-T"}


async def parse_target_discovery_query(text: str) -> dict:
    """Extract organism, disease name, and therapeutic modality from a Research Orientation Query.

    Args:
        text: Free-text research query mentioning a disease/pathogen (e.g. "말라세지아 모낭염 연고 개발")

    Returns:
        dict with keys: organism_scientific, organism_common, disease_name, modality_hint,
        context_summary (2-sentence disease background)
    """
    import json as _json
    try:
        prompt = (
            "You are a microbiologist and pharmacologist. The user is asking a research orientation "
            "query — they know a disease or pathogen but do not know which protein to target.\n\n"
            f"Query: {text}\n\n"
            "Extract the following information. If the query is in Korean, still return English values.\n"
            "Return ONLY raw JSON (no markdown, no explanation):\n"
            '{\n'
            '  "organism_scientific": "Malassezia globosa",\n'
            '  "organism_common": "Malassezia (pityrosporum) yeast",\n'
            '  "disease_name": "Malassezia folliculitis (fungal acne)",\n'
            '  "modality_hint": "topical small molecule",\n'
            '  "context_summary": "Two-sentence disease background for researchers."\n'
            '}'
        )
        raw = await _call(prompt)
        if not raw:
            return {
                "organism_scientific": "unknown",
                "organism_common": "unknown",
                "disease_name": "unknown",
                "modality_hint": "small molecule",
                "context_summary": "",
            }
        cleaned = raw.strip()
        # Strip markdown code fences if present (```json ... ``` or ``` ... ```)
        if "```" in cleaned:
            start = cleaned.find("{")
            end = cleaned.rfind("}") + 1
            if start != -1 and end > start:
                cleaned = cleaned[start:end]
        return _json.loads(cleaned)
    except Exception as exc:
        log.debug("parse_target_discovery_query failed: %s", exc)
        return {
            "organism_scientific": "unknown",
            "organism_common": "unknown",
            "disease_name": "unknown",
            "modality_hint": "small molecule",
            "context_summary": "",
        }


async def synthesize_target_discovery(
    disease_name: str,
    organism_scientific: str,
    candidates_json: str,
) -> dict:
    """Ask Gemini to synthesize top drug target candidates for a pathogen/disease.

    Returns dict with keys: ranked_rationale (str), overall_strategy (str)
    """
    import json as _json
    try:
        prompt = (
            f"You are a drug discovery expert specializing in infectious disease therapeutics.\n"
            f"Disease: {disease_name}\n"
            f"Pathogen: {organism_scientific}\n"
            f"Top protein candidates (JSON):\n{candidates_json}\n\n"
            "Provide a rigorous analysis covering:\n"
            "1. Top 3 recommended targets with mechanistic rationale\n"
            "2. Selectivity risks (human ortholog similarity warnings)\n"
            "3. Known inhibitor classes for these target families\n"
            "4. Overall therapeutic strategy (1 paragraph)\n\n"
            "Return ONLY raw JSON:\n"
            '{"ranked_rationale": "...", "overall_strategy": "..."}'
        )
        raw = await _call(prompt)
        if not raw:
            return {"ranked_rationale": "", "overall_strategy": ""}
        cleaned = raw.strip()
        if "```" in cleaned:
            start = cleaned.find("{")
            end = cleaned.rfind("}") + 1
            if start != -1 and end > start:
                cleaned = cleaned[start:end]
        return _json.loads(cleaned)
    except Exception as exc:
        log.debug("synthesize_target_discovery failed: %s", exc)
        return {"ranked_rationale": "", "overall_strategy": ""}
