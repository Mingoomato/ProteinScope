"""Synthetic Genetic Circuit Architect.

Recommends synthetic biology circuit topologies for heterologous expression or
gene-circuit integration of the target gene.  Uses iGEM BioBricks / SynBioHub
conceptual part libraries, standard regulatory logic, and Gemini 2.5 Pro for
rationale synthesis.

Design space covered
────────────────────
• Promoter tier  — constitutive (J23 series), inducible (Tet/Lac/ara/cumate)
• RBS tier       — Salis RBS Calculator predictions, Bicistronic Design (BCD)
• Logic topology — single-gene expression, Boolean gates (AND/OR/NOT), toggle
                   switches, repressilators (three-node ring oscillators)
• Terminator     — B0010/B0012 composite (iGEM), rrnB T1T2
• Vector backbone — pSB1C3 (high-copy), pSB3K3 (medium), pSB4A5 (low)
• Host chassis   — E. coli DH5α / BL21(DE3), B. subtilis 168, S. cerevisiae W303

Key citations
─────────────
• Gardner TS et al. 2000 Nature 403:339-342 doi:10.1038/35002131  (toggle switch)
• Elowitz MB & Leibler S 2000 Nature 403:335-338 doi:10.1038/35002125 (repressilator)
• Canton B et al. 2008 Nat Biotechnol 26:787-793 doi:10.1038/nbt1426 (BioBricks)
• Stanton BC et al. 2014 Science 343:1245958 doi:10.1126/science.1245958 (logic gates)
• Brophy JAN & Voigt CA 2014 Nat Methods 11:508-520 doi:10.1038/nmeth.2926 (design)
• Salis HM et al. 2009 Nat Biotechnol 27:946-950 doi:10.1038/nbt.1568 (RBS calculator)
• Lou C et al. 2012 ACS Synth Biol 1:184-190 doi:10.1021/sb200003v (bicistronic)
"""
from __future__ import annotations

import hashlib
from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

# ──────────────────────────────────────────────────────────────────────────────
# Internal part libraries (representative; real deployment would query SynBioHub)
# ──────────────────────────────────────────────────────────────────────────────

_PROMOTERS: Dict[str, Dict[str, Any]] = {
    # Constitutive — J23 Anderson series
    # Citation: Anderson 2006 iGEM Registry http://parts.igem.org/Part:BBa_J23100
    "BBa_J23100": {"type": "constitutive", "strength": "strong",   "host": "E. coli",   "sequence_hint": "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC"},
    "BBa_J23106": {"type": "constitutive", "strength": "medium",   "host": "E. coli",   "sequence_hint": "TTGACGGCTAGCTCAGTCCTAGGTATAATGCTAGC"},
    "BBa_J23117": {"type": "constitutive", "strength": "weak",     "host": "E. coli",   "sequence_hint": "TTGACAGCTAGCTCAGTCCTAGGGATTGTGCTAGC"},
    # Inducible — E. coli
    # Citation: Lutz R & Bujard H 1997 Nucleic Acids Res 25:1203 doi:10.1093/nar/25.6.1203
    "BBa_R0010": {"type": "inducible",    "inducer": "IPTG",       "host": "E. coli",   "kd_um": 1.0},
    "BBa_R0040": {"type": "inducible",    "inducer": "aTc",        "host": "E. coli",   "kd_um": 0.2},
    "BBa_K2442101": {"type": "inducible", "inducer": "arabinose",  "host": "E. coli",   "kd_um": 100.0},
    "BBa_K1319011": {"type": "inducible", "inducer": "cumate",     "host": "E. coli",   "kd_um": 0.5},
}

_RBS: Dict[str, Dict[str, Any]] = {
    # Citation: Salis HM et al. 2009 Nat Biotechnol 27:946-950 doi:10.1038/nbt.1568
    "BBa_B0034": {"strength": "strong",   "TIR_au": 10000, "note": "Reference RBS"},
    "BBa_B0032": {"strength": "medium",   "TIR_au": 3000,  "note": "Moderate translation"},
    "BBa_B0031": {"strength": "weak",     "TIR_au": 500,   "note": "Low leakage"},
    # Bicistronic Design (BCD)
    # Citation: Lou C et al. 2012 ACS Synth Biol 1:184-190 doi:10.1021/sb200003v
    "BCD2":      {"strength": "tunable",  "TIR_au": 5000,  "note": "BCD — insulates from 5′ context"},
}

_TERMINATORS: Dict[str, str] = {
    # Citation: Canton B et al. 2008 Nat Biotechnol 26:787-793 doi:10.1038/nbt1426
    "BBa_B0010": "rrnB T1 — strong rho-independent",
    "BBa_B0012": "rrnB T2 — used with B0010 for composite",
    "BBa_B0015": "B0010+B0012 composite — recommended default",
}

_BACKBONES: Dict[str, Dict[str, Any]] = {
    # Citation: Canton B et al. 2008 Nat Biotechnol 26:787-793 doi:10.1038/nbt1426
    "pSB1C3":  {"copy_number": "high",   "marker": "CmR",  "ori": "pMB1"},
    "pSB3K3":  {"copy_number": "medium", "marker": "KanR", "ori": "p15A"},
    "pSB4A5":  {"copy_number": "low",    "marker": "AmpR", "ori": "pSC101"},
}

_CIRCUIT_TOPOLOGIES: List[Dict[str, Any]] = [
    {
        "name": "Single-gene constitutive expression",
        "use_case": "stable baseline expression for structural/biochemical studies",
        "parts": ["J23100 or J23106 promoter", "BBa_B0034 RBS", "CDS", "BBa_B0015 terminator"],
        "logic": "always ON",
        "complexity": "low",
        # Citation: Brophy JAN & Voigt CA 2014 Nat Methods 11:508-520
        "reference": "Brophy 2014 doi:10.1038/nmeth.2926",
    },
    {
        "name": "Inducible T7/Lac expression (BL21-DE3)",
        "use_case": "high-yield recombinant protein production",
        "parts": ["T7lac promoter", "BBa_B0034 RBS", "CDS+His6 tag", "T7 Te terminator"],
        "logic": "IPTG-inducible",
        "complexity": "low",
        "reference": "Studier FW 2005 Protein Expr Purif 41:207 doi:10.1016/j.pep.2005.01.016",
    },
    {
        "name": "Boolean AND gate (two-input)",
        "use_case": "conditional expression requiring two simultaneous signals",
        "parts": ["Split T7 RNAP (Shis 2013)", "Input A sensor", "Input B sensor"],
        "logic": "A AND B → output",
        "complexity": "medium",
        # Citation: Stanton BC et al. 2014 Science 343:1245958
        "reference": "Stanton 2014 doi:10.1126/science.1245958",
    },
    {
        "name": "Genetic toggle switch",
        "use_case": "bistable memory; latch therapeutic gene ON/OFF state",
        "parts": ["PLtetO1 ↔ Plac/ara-1", "TetR repressor", "LacI repressor", "GOI"],
        "logic": "bistable SR-latch; IPTG → state 1; aTc → state 0",
        "complexity": "medium",
        # Citation: Gardner TS et al. 2000 Nature 403:339-342
        "reference": "Gardner 2000 doi:10.1038/35002131",
    },
    {
        "name": "Repressilator (three-node oscillator)",
        "use_case": "periodic pulsed expression; circadian mimicry",
        "parts": ["PLtetO1→lacI", "Ptrc2→cI", "PA1lacO-1→tetR", "reporter"],
        "logic": "autonomous oscillation ~150 min period in E. coli",
        "complexity": "high",
        # Citation: Elowitz MB & Leibler S 2000 Nature 403:335-338
        "reference": "Elowitz 2000 doi:10.1038/35002125",
    },
]

# Gene-size to backbone guidance
# Citation: Grieger JC & Samulski RJ 2005 J Virol 79:9933 doi:10.1128/JVI.79.15.9933
_CDS_SIZE_THRESHOLDS = {
    # CDS bp → backbone recommendation
    "<=1500": ("pSB1C3", "high-copy; fits standard BioBrick insert"),
    "<=3000": ("pSB3K3", "medium-copy; reduces metabolic burden for mid-size inserts"),
    ">3000":  ("pSB4A5", "low-copy; reduces plasmid instability for large inserts"),
}

_IGEM_REFERENCES = [
    "Gardner TS et al. 2000 Nature 403:339-342 doi:10.1038/35002131",
    "Elowitz MB & Leibler S 2000 Nature 403:335-338 doi:10.1038/35002125",
    "Canton B et al. 2008 Nat Biotechnol 26:787-793 doi:10.1038/nbt1426",
    "Stanton BC et al. 2014 Science 343:1245958 doi:10.1126/science.1245958",
    "Brophy JAN & Voigt CA 2014 Nat Methods 11:508-520 doi:10.1038/nmeth.2926",
    "Salis HM et al. 2009 Nat Biotechnol 27:946-950 doi:10.1038/nbt.1568",
    "Lou C et al. 2012 ACS Synth Biol 1:184-190 doi:10.1021/sb200003v",
]


# ──────────────────────────────────────────────────────────────────────────────
# Pydantic models
# ──────────────────────────────────────────────────────────────────────────────

class CircuitPart(BaseModel):
    """One standard biological part (BioBrick or equivalent)."""
    part_id: str
    part_type: str                          # promoter | RBS | CDS | terminator | backbone
    description: str
    key_property: Optional[str] = None     # e.g. inducer, strength, copy number


class CircuitDesignOption(BaseModel):
    """One complete circuit topology recommendation."""
    topology_name: str
    use_case: str
    logic_description: str
    complexity: str                         # low | medium | high
    parts: List[CircuitPart]
    design_notes: str
    reference: str


class SynBioHubHit(BaseModel):
    """A matching part/device found in SynBioHub/iGEM."""
    name: str
    uri: str
    description: str
    role: str                               # promoter | RBS | CDS | terminator | device


class GeneticCircuitDesign(BaseModel):
    """Full synthetic genetic circuit design for the target gene."""
    gene: str
    cds_length_bp: int
    recommended_topology: CircuitDesignOption
    alternative_topologies: List[CircuitDesignOption] = Field(default_factory=list)
    recommended_backbone: CircuitPart
    recommended_rbs: CircuitPart
    recommended_promoter: CircuitPart
    recommended_terminator: CircuitPart
    synbiohub_hits: List[SynBioHubHit] = Field(default_factory=list)
    design_rationale: str = ""
    gemini_strategy: str = ""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    references: List[str] = Field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────────────────────────────────────

def _estimate_cds_bp(sequence: str) -> int:
    """Estimate CDS length in bp from amino-acid sequence."""
    return len(sequence) * 3 + 3  # +3 for stop codon


def _choose_backbone(cds_bp: int) -> CircuitPart:
    """Select backbone by CDS size.

    Citation: Grieger JC & Samulski RJ 2005 J Virol 79:9933 doi:10.1128/JVI.79.15.9933
    (plasmid copy number → metabolic burden relationship)
    """
    if cds_bp <= 1500:
        key, note = "pSB1C3", "high-copy; standard BioBrick chassis"
    elif cds_bp <= 3000:
        key, note = "pSB3K3", "medium-copy; lower burden for mid-size inserts"
    else:
        key, note = "pSB4A5", "low-copy; minimises instability for large CDS"
    info = _BACKBONES[key]
    return CircuitPart(
        part_id=key,
        part_type="backbone",
        description=f"{info['ori']} ori, {info['marker']}, {info['copy_number']}-copy",
        key_property=note,
    )


def _choose_promoter(function_description: str, cds_bp: int) -> CircuitPart:
    """Select promoter based on inferred expression strategy.

    Heuristics (Brophy 2014 doi:10.1038/nmeth.2926):
    • Toxic / tightly regulated proteins → inducible (IPTG or aTc)
    • Large CDS (> 2000 bp) → medium-strength constitutive to reduce burden
    • Default → strong constitutive J23100
    """
    desc_lower = (function_description or "").lower()
    # inducibility heuristic — toxicity or regulation key-words
    toxic_keywords = ["kinase", "protease", "toxin", "transcription factor", "polymerase"]
    is_regulated = any(kw in desc_lower for kw in toxic_keywords)

    if is_regulated:
        pid, info = "BBa_R0040", _PROMOTERS["BBa_R0040"]
        desc = f"aTc-inducible (TetR-repressed); Kd ≈ {info['kd_um']} µM aTc"
        prop = "Recommended for toxic/tightly-regulated proteins — inducible"
    elif cds_bp > 2000:
        pid, info = "BBa_J23106", _PROMOTERS["BBa_J23106"]
        desc = "Constitutive medium-strength; reduces metabolic burden for large CDS"
        prop = "Medium constitutive (J23106 Anderson series)"
    else:
        pid, info = "BBa_J23100", _PROMOTERS["BBa_J23100"]
        desc = "Constitutive strong; reference Anderson promoter"
        prop = "Strong constitutive (J23100 Anderson series)"

    return CircuitPart(part_id=pid, part_type="promoter", description=desc, key_property=prop)


def _choose_rbs(cds_bp: int) -> CircuitPart:
    """Select RBS.

    Citation: Salis HM et al. 2009 Nat Biotechnol 27:946-950 doi:10.1038/nbt.1568
    • Large CDS (> 3000 bp) → BCD2 for context-insulation
    • Else → BBa_B0034 (reference RBS, TIR 10000 AU)
    """
    if cds_bp > 3000:
        return CircuitPart(
            part_id="BCD2",
            part_type="RBS",
            description="Bicistronic Design 2 — insulates RBS strength from 5′ context",
            key_property="TIR ~5000 AU (Lou 2012 doi:10.1021/sb200003v)",
        )
    return CircuitPart(
        part_id="BBa_B0034",
        part_type="RBS",
        description="Reference RBS; TIR ~10 000 AU (Salis 2009 doi:10.1038/nbt.1568)",
        key_property="Strong translation initiation",
    )


def _choose_terminator() -> CircuitPart:
    """Always recommend composite B0015.

    Citation: Canton B et al. 2008 Nat Biotechnol 26:787-793 doi:10.1038/nbt1426
    """
    return CircuitPart(
        part_id="BBa_B0015",
        part_type="terminator",
        description="B0010+B0012 composite rho-independent terminator — recommended default",
        key_property="Termination efficiency >99% in E. coli",
    )


def _build_topology(
    topo_dict: Dict[str, Any],
    promoter: CircuitPart,
    rbs: CircuitPart,
    terminator: CircuitPart,
    gene: str,
) -> CircuitDesignOption:
    """Build a CircuitDesignOption from a topology template + selected parts."""
    parts = [promoter, rbs]
    parts.append(CircuitPart(part_id="CDS", part_type="CDS", description=f"{gene} coding sequence"))
    parts.append(terminator)
    return CircuitDesignOption(
        topology_name=topo_dict["name"],
        use_case=topo_dict["use_case"],
        logic_description=topo_dict["logic"],
        complexity=topo_dict["complexity"],
        parts=parts,
        design_notes=" | ".join(topo_dict.get("parts", [])),
        reference=topo_dict["reference"],
    )


def _seed_topology_choice(gene: str) -> int:
    """Deterministic secondary topology choice seeded by gene name hash."""
    # Ensures same gene always gets same alternative topology (reproducibility)
    h = int(hashlib.sha256(gene.encode()).hexdigest(), 16)
    return h % len(_CIRCUIT_TOPOLOGIES)


async def _mock_synbiohub_lookup(gene: str) -> List[SynBioHubHit]:
    """Mock SynBioHub lookup — returns conceptual hits.

    Real deployment would GET:
    https://synbiohub.org/search?q={gene}&objectType=ComponentDefinition
    Citation: McLaughlin JA et al. 2018 ACS Synth Biol 7:2745 doi:10.1021/acssynbio.8b00140
    """
    # Deterministic stub — in production replace with actual SynBioHub REST call
    return [
        SynBioHubHit(
            name=f"{gene}_expression_cassette (conceptual)",
            uri=f"https://synbiohub.org/search?q={gene}",
            description=f"Search SynBioHub for existing {gene} expression devices",
            role="device",
        ),
        SynBioHubHit(
            name=f"iGEM Registry: {gene} CDS parts",
            uri=f"http://parts.igem.org/cgi/search.cgi?q={gene}",
            description=f"iGEM Registry search for {gene}-related BioBricks",
            role="CDS",
        ),
    ]


# ──────────────────────────────────────────────────────────────────────────────
# Public entry-point
# ──────────────────────────────────────────────────────────────────────────────

async def run_genetic_circuit_design(
    gene: str,
    sequence: str,
    function_description: str = "",
    step_cb=None,
) -> Optional[GeneticCircuitDesign]:
    """Design a synthetic genetic circuit for the target gene.

    Parameters
    ----------
    gene:
        HGNC gene symbol.
    sequence:
        Canonical amino-acid sequence (single-letter code).
    function_description:
        Plain-text UniProt function annotation (guides promoter choice heuristics).
    step_cb:
        Optional async callable for SSE progress reporting.

    Returns
    -------
    GeneticCircuitDesign or None on any unrecoverable error.
    """
    try:
        # ── Step 1: estimate CDS size ─────────────────────────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "estimating_cds_size", 0.1)

        cds_bp = _estimate_cds_bp(sequence) if sequence else 750  # fallback ~250 aa

        # ── Step 2: select parts ──────────────────────────────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "selecting_parts", 0.25)

        backbone   = _choose_backbone(cds_bp)
        promoter   = _choose_promoter(function_description, cds_bp)
        rbs        = _choose_rbs(cds_bp)
        terminator = _choose_terminator()

        # ── Step 3: primary topology — single-gene constitutive / inducible ───
        if step_cb:
            await step_cb("genetic_circuit", "selecting_topology", 0.45)

        # Primary topology: choose by promoter type
        if "inducible" in promoter.key_property.lower():
            primary_topo = _CIRCUIT_TOPOLOGIES[1]  # IPTG/aTc inducible T7
        else:
            primary_topo = _CIRCUIT_TOPOLOGIES[0]  # constitutive

        recommended = _build_topology(primary_topo, promoter, rbs, terminator, gene)

        # ── Step 4: alternative topologies ───────────────────────────────────
        alt_idx = _seed_topology_choice(gene)
        alt_topos: List[CircuitDesignOption] = []
        for i, t in enumerate(_CIRCUIT_TOPOLOGIES[2:], start=2):  # skip basic ones
            if i == alt_idx % len(_CIRCUIT_TOPOLOGIES) or len(alt_topos) < 2:
                alt_topos.append(_build_topology(t, promoter, rbs, terminator, gene))
            if len(alt_topos) >= 2:
                break

        # ── Step 5: SynBioHub / iGEM conceptual lookup ───────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "querying_synbiohub", 0.6)

        synbiohub_hits = await _mock_synbiohub_lookup(gene)

        # ── Step 6: design rationale text ────────────────────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "writing_rationale", 0.75)

        rationale_parts = [
            f"CDS estimated at {cds_bp} bp ({len(sequence) if sequence else '~250'} aa).",
            f"Backbone: {backbone.part_id} ({backbone.key_property}).",
            f"Promoter: {promoter.part_id} — {promoter.key_property}.",
            f"RBS: {rbs.part_id} — {rbs.key_property}.",
            f"Terminator: {terminator.part_id} — {terminator.description}.",
            f"Primary topology: {recommended.topology_name} ({recommended.complexity} complexity).",
        ]
        design_rationale = "  ".join(rationale_parts)

        # ── Step 7: Gemini synthesis ──────────────────────────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "gemini_synthesis", 0.88)

        gemini_strategy = ""
        try:
            from core.gemini_interpreter import _call  # type: ignore
            prompt = (
                f"You are a synthetic biology expert.  Design a genetic circuit for the gene "
                f"'{gene}' ({len(sequence) if sequence else 'unknown'} aa).\n\n"
                f"Protein function: {function_description or 'not provided'}\n\n"
                f"Selected parts:\n"
                f"  Promoter  : {promoter.part_id} — {promoter.description}\n"
                f"  RBS       : {rbs.part_id} — {rbs.description}\n"
                f"  Backbone  : {backbone.part_id} — {backbone.description}\n"
                f"  Topology  : {recommended.topology_name}\n\n"
                f"Please provide:\n"
                f"1. Why this circuit topology is appropriate for {gene}\n"
                f"2. Potential metabolic-burden or leaky-expression concerns\n"
                f"3. Recommended chassis organism and culture conditions\n"
                f"4. Alternative circuit strategies if primary fails\n"
                f"5. Key validation experiments (Western blot, qRT-PCR, flow cytometry)\n\n"
                f"Write at the level of a synthetic biology PhD student."
            )
            gemini_strategy = await _call(prompt)
        except Exception:
            gemini_strategy = ""

        # ── Step 8: assemble result ───────────────────────────────────────────
        if step_cb:
            await step_cb("genetic_circuit", "complete", 1.0)

        return GeneticCircuitDesign(
            gene=gene.upper(),
            cds_length_bp=cds_bp,
            recommended_topology=recommended,
            alternative_topologies=alt_topos,
            recommended_backbone=backbone,
            recommended_rbs=rbs,
            recommended_promoter=promoter,
            recommended_terminator=terminator,
            synbiohub_hits=synbiohub_hits,
            design_rationale=design_rationale,
            gemini_strategy=gemini_strategy or "",
            timestamp=datetime.utcnow(),
            references=_IGEM_REFERENCES,
        )

    except Exception:
        return None
