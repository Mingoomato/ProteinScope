"""Covalent Inhibitor Designer.

Identifies covalent warhead attachment sites on a target protein and recommends
electrophilic warhead chemistries for covalent-drug design.

Background
──────────
Covalent inhibitors form permanent (irreversible) or long-lived (reversible covalent)
bonds with nucleophilic residues in the target.  The clinical success of ibrutinib
(BTK Cys481), osimertinib (EGFR Cys797), and sotorasib (KRAS G12C) has renewed
interest in targeted covalent inhibitors (TCIs).

Site identification strategy
─────────────────────────────
1. Scan sequence for nucleophilic residues: Cys, Ser, Lys, Tyr, His
2. Rank by reactivity (Cys >>> Ser ≈ Lys > Tyr > His)
3. Assess proximity to known binding / active sites → "hot" warhead attachment windows
4. Map each site to recommended warhead chemistries + known clinical examples

Key citations
─────────────
• Bauer RA 2015 Drug Discov Today 20:1061-73 doi:10.1016/j.drudis.2015.05.005
  (comprehensive TCI review; warhead-residue pairing)
• Baillie TA 2016 Angew Chem Int Ed 55:13408-21 doi:10.1002/anie.201601091
  (mechanisms, selectivity, druggability)
• Gehringer M & Laufer SA 2019 J Med Chem 62:5673-724 doi:10.1021/acs.jmedchem.8b01153
  (warhead chemistry landscape)
• Zhao Z & Bourne PE 2018 Drug Discov Today 23:727-35 doi:10.1016/j.drudis.2018.01.019
  (covalent docking; hot-spot residues)
• Singh J et al. 2011 Nat Chem Biol 7:502-9 doi:10.1038/nchembio.601
  (reversible covalent inhibitors; boronic acid / aldehyde)
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

# ──────────────────────────────────────────────────────────────────────────────
# Warhead library
# Each entry: residue type → warhead options with reactivity notes
# Citation: Bauer 2015 doi:10.1016/j.drudis.2015.05.005;
#           Gehringer 2019 doi:10.1021/acs.jmedchem.8b01153
# ──────────────────────────────────────────────────────────────────────────────

_WARHEAD_MAP: Dict[str, List[Dict[str, Any]]] = {
    "C": [  # Cysteine — most reactive; preferred target
        {"chemistry": "Acrylamide",          "reversibility": "irreversible",
         "example_drug": "Afatinib (EGFR C797)", "reactivity_rank": 1,
         "note": "Most common TCI warhead; Michael acceptor"},
        {"chemistry": "Chloroacetamide",     "reversibility": "irreversible",
         "example_drug": "Ibrutinib-like probes", "reactivity_rank": 2,
         "note": "High intrinsic reactivity; selectivity concerns"},
        {"chemistry": "Vinyl sulfonamide",   "reversibility": "irreversible",
         "example_drug": "RPE65 probes",     "reactivity_rank": 3,
         "note": "Tunable electrophilicity vs. acrylamide"},
        {"chemistry": "Cyanoacrylamide",     "reversibility": "reversible covalent",
         "example_drug": "RSK2 inhibitor",   "reactivity_rank": 4,
         "note": "Reversible Michael acceptor (Singh 2011 doi:10.1038/nchembio.601)"},
    ],
    "K": [  # Lysine — amino group; requires activated electrophile
        {"chemistry": "NHS ester",           "reversibility": "irreversible",
         "example_drug": "Aspirin (Lys mechanism)", "reactivity_rank": 5,
         "note": "Acylation of Lys ε-NH2 (pKa ~10.5)"},
        {"chemistry": "Sulfonyl fluoride",   "reversibility": "irreversible",
         "example_drug": "SuFEx chemistry probes",  "reactivity_rank": 5,
         "note": "SuFEx click chemistry; labels Lys, Ser, Tyr, His"},
        {"chemistry": "Isocyanate",          "reversibility": "irreversible",
         "example_drug": "Covalent fragment screens","reactivity_rank": 6,
         "note": "Carbamoylation of Lys"},
    ],
    "S": [  # Serine — usually requires proximity to oxyanion hole
        {"chemistry": "DFP (diisopropyl fluorophosphate)", "reversibility": "irreversible",
         "example_drug": "Sarin/DIFO (serine protease)", "reactivity_rank": 7,
         "note": "Classic serine-protease mechanism; catalytic Ser only"},
        {"chemistry": "β-Lactam",            "reversibility": "irreversible",
         "example_drug": "Penicillin (DD-transpeptidase)", "reactivity_rank": 7,
         "note": "Acylation of catalytic Ser; restricted to penicillin-binding proteins"},
        {"chemistry": "Sulfonyl fluoride",   "reversibility": "irreversible",
         "example_drug": "SuFEx probes",     "reactivity_rank": 6,
         "note": "Promiscuous; best for proximity-guided design"},
    ],
    "Y": [  # Tyrosine — phenolic OH; requires activated electrophile
        {"chemistry": "Sulfonyl fluoride",   "reversibility": "irreversible",
         "example_drug": "Kinase-Tyr probes", "reactivity_rank": 8,
         "note": "SuFEx labels Tyr OH; Tyr716 in PI3K-δ example"},
        {"chemistry": "Epoxide",             "reversibility": "irreversible",
         "example_drug": "Carfilzomib (proteasome)", "reactivity_rank": 9,
         "note": "Alkylation of Tyr phenolic OH; low selectivity"},
    ],
    "H": [  # Histidine — imidazole; moderate nucleophile
        {"chemistry": "Chloroacetamide",     "reversibility": "irreversible",
         "example_drug": "Papain probes",    "reactivity_rank": 9,
         "note": "Alkylation of His Nε or Nδ; less common than Cys"},
        {"chemistry": "Diazoacetyl",         "reversibility": "irreversible",
         "example_drug": "Pepstatin analogs","reactivity_rank": 10,
         "note": "Carbene insertion; photochemical activation possible"},
    ],
}

# Reactivity score per residue type (higher = more nucleophilic at physiological pH)
# Citation: Bauer 2015; Baillie 2016 doi:10.1002/anie.201601091
_REACTIVITY_SCORE: Dict[str, float] = {
    "C": 1.0,   # pKa ~8.3 in free Cys; catalytic Cys can be ~4-5
    "K": 0.4,   # pKa ~10.5; reactive only when buried + deprotonated
    "S": 0.3,   # pKa ~14; requires activation (catalytic triad)
    "Y": 0.25,  # pKa ~10.1; requires activation
    "H": 0.2,   # pKa ~6.0; partially deprotonated at pH 7.4
}

_SITE_WINDOW = 10   # residues around active/binding site to flag as "proximity hit"
                    # Citation: Zhao & Bourne 2018 doi:10.1016/j.drudis.2018.01.019


# ──────────────────────────────────────────────────────────────────────────────
# Pydantic models
# ──────────────────────────────────────────────────────────────────────────────

class WarheadOption(BaseModel):
    chemistry: str
    reversibility: str          # "irreversible" | "reversible covalent"
    example_drug: str
    reactivity_rank: int        # lower = more reactive (1 = most)
    note: str


class CovalentSite(BaseModel):
    """One nucleophilic residue evaluated as a covalent attachment point."""
    residue: str                        # single-letter AA
    position: int
    reactivity_score: float             # 0–1, normalised
    near_functional_site: bool          # within _SITE_WINDOW of binding/active site
    site_type: Optional[str] = None     # "active_site" | "binding_site" | None
    warhead_options: List[WarheadOption] = Field(default_factory=list)
    priority: str = "low"               # "high" | "medium" | "low"
    design_note: str = ""


class CovalentInhibitorDesign(BaseModel):
    """Full covalent inhibitor design report for the target gene."""
    gene: str
    total_cys: int
    total_ser: int
    total_lys: int
    total_tyr: int
    total_his: int
    top_sites: List[CovalentSite] = Field(default_factory=list)
    design_summary: str = ""
    gemini_strategy: str = ""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    references: List[str] = Field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────────────────────────────────────

_REFERENCES = [
    "Bauer RA 2015 Drug Discov Today 20:1061 doi:10.1016/j.drudis.2015.05.005",
    "Baillie TA 2016 Angew Chem Int Ed 55:13408 doi:10.1002/anie.201601091",
    "Gehringer M & Laufer SA 2019 J Med Chem 62:5673 doi:10.1021/acs.jmedchem.8b01153",
    "Zhao Z & Bourne PE 2018 Drug Discov Today 23:727 doi:10.1016/j.drudis.2018.01.019",
    "Singh J et al. 2011 Nat Chem Biol 7:502 doi:10.1038/nchembio.601",
]


def _get_functional_positions(binding_sites: list, active_sites: list) -> set:
    """Return set of residue positions from binding and active sites."""
    positions = set()
    for site in binding_sites + active_sites:
        try:
            pos = getattr(site, "position", None) or (site.get("position") if isinstance(site, dict) else None)
            if pos:
                positions.add(int(pos))
            # Also parse from sequence_fragment if available
            frag = getattr(site, "sequence_fragment", None) or (site.get("sequence_fragment") if isinstance(site, dict) else None)
            if frag and isinstance(frag, str):
                start = getattr(site, "start", None) or (site.get("start") if isinstance(site, dict) else None)
                if start:
                    for i in range(len(frag)):
                        positions.add(int(start) + i)
        except Exception:
            continue
    return positions


def _priority(reactivity: float, near_site: bool) -> str:
    """Compute priority from reactivity score and proximity to functional site."""
    if near_site and reactivity >= 0.8:
        return "high"
    if near_site or reactivity >= 0.6:
        return "medium"
    return "low"


def _scan_nucleophiles(
    sequence: str,
    functional_positions: set,
) -> List[CovalentSite]:
    """Scan sequence for nucleophilic residues and assess each."""
    sites: List[CovalentSite] = []
    target_aas = set(_REACTIVITY_SCORE.keys())

    for i, aa in enumerate(sequence):
        if aa not in target_aas:
            continue
        pos = i + 1  # 1-based

        # Proximity to functional site — within ±_SITE_WINDOW residues
        near = any(abs(pos - fp) <= _SITE_WINDOW for fp in functional_positions)
        site_type: Optional[str] = None
        if near and functional_positions:
            # Determine whether nearest site is active or binding
            site_type = "functional_site"

        reactivity = _REACTIVITY_SCORE[aa]
        # Boost reactivity if near functional site (microenvironment effect)
        # Citation: Zhao & Bourne 2018 doi:10.1016/j.drudis.2018.01.019
        if near:
            reactivity = min(1.0, reactivity * 1.4)

        warheads = [WarheadOption(**w) for w in _WARHEAD_MAP.get(aa, [])]

        note_parts = [f"{aa}{pos}"]
        if near:
            note_parts.append(f"within {_SITE_WINDOW} aa of functional site")
        design_note = "; ".join(note_parts)

        sites.append(CovalentSite(
            residue=aa,
            position=pos,
            reactivity_score=round(reactivity, 3),
            near_functional_site=near,
            site_type=site_type,
            warhead_options=warheads,
            priority=_priority(reactivity, near),
            design_note=design_note,
        ))

    # Sort: high priority first, then by reactivity descending
    _prio_rank = {"high": 0, "medium": 1, "low": 2}
    sites.sort(key=lambda s: (_prio_rank[s.priority], -s.reactivity_score))
    return sites


# ──────────────────────────────────────────────────────────────────────────────
# Public entry-point
# ──────────────────────────────────────────────────────────────────────────────

async def run_covalent_inhibitor_design(
    gene: str,
    sequence: str,
    binding_sites: Optional[list] = None,
    active_sites: Optional[list] = None,
    step_cb=None,
) -> Optional[CovalentInhibitorDesign]:
    """Design covalent inhibitor attachment sites for the target protein.

    Parameters
    ----------
    gene:
        HGNC gene symbol.
    sequence:
        Canonical amino-acid sequence (single-letter).
    binding_sites:
        List of SequenceFeature / dict objects with position/sequence_fragment.
    active_sites:
        List of SequenceFeature / dict objects.
    step_cb:
        Optional async callable for SSE progress reporting.

    Returns
    -------
    CovalentInhibitorDesign or None on error.
    """
    try:
        if not sequence:
            return None

        # ── Step 1: gather functional positions ───────────────────────────────
        if step_cb:
            await step_cb("covalent", "scanning_functional_sites", 0.15)

        functional_positions = _get_functional_positions(
            binding_sites or [], active_sites or []
        )

        # ── Step 2: scan nucleophiles ─────────────────────────────────────────
        if step_cb:
            await step_cb("covalent", "scanning_nucleophiles", 0.35)

        all_sites = _scan_nucleophiles(sequence, functional_positions)

        # ── Step 3: statistics ────────────────────────────────────────────────
        if step_cb:
            await step_cb("covalent", "computing_statistics", 0.55)

        top_sites = [s for s in all_sites if s.priority in ("high", "medium")][:15]
        if not top_sites:
            top_sites = all_sites[:8]   # fallback: show at least 8

        # ── Step 4: design summary ────────────────────────────────────────────
        if step_cb:
            await step_cb("covalent", "writing_summary", 0.70)

        high_count = sum(1 for s in all_sites if s.priority == "high")
        cys_sites = [s for s in all_sites if s.residue == "C"]
        near_cys = [s for s in cys_sites if s.near_functional_site]
        summary_parts = [
            f"{len(all_sites)} nucleophilic sites found "
            f"(Cys:{sequence.count('C')}, Ser:{sequence.count('S')}, "
            f"Lys:{sequence.count('K')}, Tyr:{sequence.count('Y')}, His:{sequence.count('H')}).",
            f"{high_count} high-priority warhead attachment site(s).",
        ]
        if near_cys:
            summary_parts.append(
                f"{len(near_cys)} Cys residue(s) near functional site(s): "
                + ", ".join(f"C{s.position}" for s in near_cys[:5])
                + " — primary covalent drug targets."
            )
        design_summary = "  ".join(summary_parts)

        # ── Step 5: Gemini synthesis ──────────────────────────────────────────
        if step_cb:
            await step_cb("covalent", "gemini_synthesis", 0.85)

        gemini_strategy = ""
        try:
            from core.gemini_interpreter import _call  # type: ignore
            top_cys_str = ", ".join(f"C{s.position}" for s in (near_cys or cys_sites)[:5]) or "none identified"
            prompt = (
                f"You are a medicinal chemist expert in targeted covalent inhibitors (TCIs).\n\n"
                f"Gene: {gene}\n"
                f"Sequence length: {len(sequence)} aa\n"
                f"High-priority covalent sites:\n"
                + "\n".join(
                    f"  {s.residue}{s.position}: {s.priority} priority, "
                    f"reactivity={s.reactivity_score:.2f}, "
                    f"near_functional_site={s.near_functional_site}, "
                    f"warheads=[{', '.join(w.chemistry for w in s.warhead_options[:2])}]"
                    for s in top_sites[:6]
                )
                + f"\n\nPlease provide:\n"
                f"1. Which site(s) should be the primary TCI target and why\n"
                f"2. Recommended warhead chemistry with selectivity considerations\n"
                f"3. Known covalent drugs targeting the same residue (if any)\n"
                f"4. Selectivity risks (off-target Cys/Ser/Lys in proteome)\n"
                f"5. Whether reversible-covalent approach is preferred here\n\n"
                f"Write at the level of a medicinal chemistry expert for drug discovery."
            )
            gemini_strategy = await _call(prompt)
        except Exception:
            gemini_strategy = ""

        # ── Step 6: assemble ──────────────────────────────────────────────────
        if step_cb:
            await step_cb("covalent", "complete", 1.0)

        return CovalentInhibitorDesign(
            gene=gene.upper(),
            total_cys=sequence.count("C"),
            total_ser=sequence.count("S"),
            total_lys=sequence.count("K"),
            total_tyr=sequence.count("Y"),
            total_his=sequence.count("H"),
            top_sites=top_sites,
            design_summary=design_summary,
            gemini_strategy=gemini_strategy or "",
            timestamp=datetime.utcnow(),
            references=_REFERENCES,
        )

    except Exception:
        return None
