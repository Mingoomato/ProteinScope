"""Metalloenzyme Analyzer — Metal coordination geometry + PROPKA pKa prediction.

Background — Metalloenzymes:
─────────────────────────────────────────────────────────────────────────────
Metalloenzymes are enzymes that require metal ion cofactors for catalytic
activity.  The metal performs one or more roles:

  Lewis acid catalysis  — the metal activates substrates by polarising bonds
                          (e.g. Zn²⁺ in carbonic anhydrase, carboxypeptidase)
  Redox cycling         — one-electron transfers via Fe²⁺/Fe³⁺ or Cu⁺/Cu²⁺
                          (e.g. cytochrome P450, laccase)
  Structural scaffold   — stabilising active-site geometry (e.g. Zn²⁺ in
                          zinc-finger domains)
  Electrophile activation — Mg²⁺ or Mn²⁺ in nucleotide-processing enzymes

Metal coordination geometry:
  Tetrahedral  (4-coord)  — common for Zn²⁺ (His₂CysGlu, His₃water, Cys₄)
  Square planar (4-coord)  — Cu²⁺, Ni²⁺
  Octahedral   (6-coord)  — Fe²⁺/Fe³⁺, Mg²⁺, Ca²⁺, Mn²⁺
  Trigonal bipyramidal    — intermediate during catalysis

pKa perturbation (PROPKA):
  The electrostatic environment near a metal centre dramatically shifts pKa
  values of catalytic residues.  PROPKA3 (Olsson et al. 2011) calculates
  residue pKa values from a PDB structure.  Shifts > 1 pH unit from the
  model-compound pKa are considered functionally significant.

  Model-compound pKa reference values:
    HIS  6.5   CYS  9.0   ASP  3.8   GLU  4.1
    LYS 10.5   ARG 12.5   TYR 10.5

Inhibition strategies for metalloenzymes:
  1. Metal chelators     — competitive or irreversible removal of the metal ion
                           (e.g. EDTA, hydroxamates for Zn MMPs, deferiprone for Fe)
  2. Competitive inhibitors — mimic transition state at the active site
  3. Allosteric inhibitors  — bind distant site, alter coordination geometry

This analyzer:
  1. Detects metal ions in the protein's cofactor list
  2. Parses UniProt binding-site features to identify coordinating residues
  3. Downloads the AlphaFold PDB structure
  4. Runs PROPKA3 to predict pKa shifts of active-site residues
  5. Calls Gemini to synthesise the catalytic mechanism + inhibition strategies
  6. Assembles a MetalloenzymeAnalysis report
"""

from __future__ import annotations

import asyncio
import json as _json
import re
import tempfile
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade, DataProvenance  # noqa: F401


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class MetalCoordination(BaseModel):
    metal_ion: str                                 # e.g. "Zn2+"
    coordinating_residues: List[str]               # e.g. ["H94", "H96", "H119", "E117"]
    coordination_number: int = 0
    geometry: str = ""                             # "tetrahedral", "octahedral", etc.
    source: str = "UniProt"


class PKaShift(BaseModel):
    residue: str                                   # e.g. "H119"
    position: int
    model_compound_pka: float                      # typical pKa for this AA type
    predicted_pka: float                           # PROPKA prediction
    shift: float                                   # predicted - model_compound
    functional_significance: str = ""             # "catalytic proton shuttle", etc.


class MetalloenzymeAnalysis(BaseModel):
    gene: str
    metal_ions: List[str] = Field(default_factory=list)           # e.g. ["Zn2+"]
    coordination_sites: List[MetalCoordination] = Field(default_factory=list)
    pka_shifts: List[PKaShift] = Field(default_factory=list)
    propka_available: bool = False
    catalytic_mechanism: str = ""
    inhibition_strategies: List[str] = Field(default_factory=list)
    gemini_synthesis: str = ""
    provenance: Optional[DataProvenance] = None


# ---------------------------------------------------------------------------
# Metal detection helpers
# ---------------------------------------------------------------------------

# Canonical metal names → canonical ion label
_METAL_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"\bzinc\b",       re.I), "Zn2+"),
    (re.compile(r"\bfe\b|iron",    re.I), "Fe2+/Fe3+"),
    (re.compile(r"\bcopper\b",     re.I), "Cu2+"),
    (re.compile(r"\bmanganese\b",  re.I), "Mn2+"),
    (re.compile(r"\bmagnesium\b",  re.I), "Mg2+"),
    (re.compile(r"\bcalcium\b",    re.I), "Ca2+"),
    (re.compile(r"\bcobalt\b",     re.I), "Co2+"),
    (re.compile(r"\bnickel\b",     re.I), "Ni2+"),
    (re.compile(r"\bmolybdenum\b", re.I), "Mo"),
]

def _detect_metals(cofactors: list[str]) -> list[str]:
    """Return list of canonical ion labels found in the cofactor list."""
    found: list[str] = []
    seen: set[str] = set()
    combined = " ".join(cofactors)
    for pattern, label in _METAL_PATTERNS:
        if pattern.search(combined) and label not in seen:
            found.append(label)
            seen.add(label)
    return found


# ---------------------------------------------------------------------------
# UniProt binding-site parser
# ---------------------------------------------------------------------------

def _parse_uniprot_binding_sites(entry: dict, metal_ions: list[str]) -> list[MetalCoordination]:
    """Extract metal-coordinating residues from UniProt feature annotations."""
    if not isinstance(entry, dict):
        return []

    # Build a set of metal keywords to match against ligandDescription
    metal_keywords = set()
    for ion in metal_ions:
        # e.g. "Zn2+" → "zn", "Fe2+/Fe3+" → "fe", "iron"
        stem = re.sub(r"[^a-zA-Z]", "", ion.split("/")[0]).lower()
        metal_keywords.add(stem)

    # Also map ion labels back to plain names for matching
    _extra_keywords = {
        "zn": ["zinc", "zn"],
        "fe": ["iron", "fe", "heme", "haem"],
        "cu": ["copper", "cu"],
        "mn": ["manganese", "mn"],
        "mg": ["magnesium", "mg"],
        "ca": ["calcium", "ca"],
        "co": ["cobalt", "co"],
        "ni": ["nickel", "ni"],
        "mo": ["molybdenum", "mo"],
    }
    expanded: set[str] = set()
    for stem in metal_keywords:
        expanded.update(_extra_keywords.get(stem, [stem]))
    metal_keywords = expanded

    # Group residues by ligandDescription (= one site per distinct metal ligand)
    site_map: dict[str, list[str]] = {}
    for feat in entry.get("features", []):
        if feat.get("type") != "Binding site":
            continue
        ligand = feat.get("ligandDescription", feat.get("ligand", {}).get("name", "")) or ""
        ligand_lower = ligand.lower()
        if not any(kw in ligand_lower for kw in metal_keywords):
            continue
        # Derive position
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        pos = start if start is not None else end
        # Derive one-letter or three-letter residue code
        seq_desc = feat.get("description", "")
        # e.g. description may say "Zinc; via tele nitrogen"
        aa_code = ""
        # Try to get from featureId or description
        # Fall back to position label
        if pos is not None:
            aa_code = str(pos)
        if not aa_code:
            continue
        ligand_key = ligand or "metal"
        site_map.setdefault(ligand_key, []).append(aa_code)

    if not site_map:
        return []

    # Convert to MetalCoordination objects
    coordinations: list[MetalCoordination] = []
    for ligand_key, residues in site_map.items():
        cn = len(residues)
        # Infer geometry from coordination number
        if cn == 4:
            geometry = "tetrahedral"
        elif cn == 5:
            geometry = "trigonal bipyramidal"
        elif cn == 6:
            geometry = "octahedral"
        else:
            geometry = f"{cn}-coordinate"

        # Determine which metal ion this site belongs to
        ion_label = "metal"
        for ion in metal_ions:
            stem = re.sub(r"[^a-zA-Z]", "", ion.split("/")[0]).lower()
            kws = _extra_keywords.get(stem, [stem])
            if any(kw in ligand_key.lower() for kw in kws):
                ion_label = ion
                break

        coordinations.append(MetalCoordination(
            metal_ion=ion_label,
            coordinating_residues=residues,
            coordination_number=cn,
            geometry=geometry,
            source="UniProt",
        ))

    return coordinations


# ---------------------------------------------------------------------------
# PROPKA output parser
# ---------------------------------------------------------------------------

# Model-compound pKa values (from Thurlkill et al. 2006 / PROPKA defaults)
_MODEL_PKA: dict[str, float] = {
    "HIS": 6.5,
    "CYS": 9.0,
    "ASP": 3.8,
    "GLU": 4.1,
    "LYS": 10.5,
    "ARG": 12.5,
    "TYR": 10.5,
}

# PROPKA output line pattern:
#  HIS  94 A    6.43      X   3 (1.0)     MET  96 A     0.0  (    -)
_PROPKA_LINE_RE = re.compile(
    r"^\s*([A-Z]{3})\s+(\d+)\s+\S+\s+([\d.]+)\s"
)

def _parse_propka_output(pka_text: str, shift_threshold: float = 1.0) -> list[PKaShift]:
    """Parse PROPKA .pka output and return residues with significant pKa shifts."""
    shifts: list[PKaShift] = []
    in_summary = False

    for line in pka_text.splitlines():
        # The summary table starts after "SUMMARY OF THIS PREDICTION"
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        # Stop at separator lines or blank sections
        if line.startswith("---") or line.strip() == "":
            continue

        m = _PROPKA_LINE_RE.match(line)
        if not m:
            continue

        res_type = m.group(1).upper()
        res_num = int(m.group(2))
        predicted_pka = float(m.group(3))

        model_pka = _MODEL_PKA.get(res_type)
        if model_pka is None:
            continue  # Only track titratable residues with known model pKa

        shift = predicted_pka - model_pka
        if abs(shift) < shift_threshold:
            continue

        # Annotate functional significance for large shifts
        significance = ""
        if res_type == "HIS" and predicted_pka > 7.0:
            significance = "Potential catalytic proton shuttle (upshifted His)"
        elif res_type == "HIS" and predicted_pka < 5.5:
            significance = "Metal-coordinated His with suppressed pKa"
        elif res_type == "CYS" and predicted_pka < 7.0:
            significance = "Thiolate nucleophile — metal-activated Cys"
        elif res_type == "ASP" and predicted_pka > 6.0:
            significance = "Buried/charge-neutralised Asp — possible proton relay"
        elif res_type == "GLU" and predicted_pka > 6.5:
            significance = "Buried Glu — likely electrostatic stabilisation role"
        elif res_type == "LYS" and predicted_pka < 8.5:
            significance = "Downshifted Lys — potential catalytic base"

        shifts.append(PKaShift(
            residue=f"{res_type[0]}{res_num}",
            position=res_num,
            model_compound_pka=model_pka,
            predicted_pka=predicted_pka,
            shift=round(shift, 2),
            functional_significance=significance,
        ))

    return shifts


# ---------------------------------------------------------------------------
# PROPKA runner
# ---------------------------------------------------------------------------

async def _run_propka(pdb_path: str) -> tuple[bool, list[PKaShift]]:
    """Run PROPKA3 on the given PDB file.  Returns (available, pka_shifts)."""
    try:
        proc = await asyncio.create_subprocess_exec(
            "propka3", pdb_path,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=120)

        # PROPKA writes the .pka file next to the input PDB
        import os
        pka_path = os.path.splitext(pdb_path)[0] + ".pka"
        if not os.path.exists(pka_path):
            # Some versions write to stdout
            pka_text = stdout.decode("utf-8", errors="replace")
        else:
            with open(pka_path, "r", encoding="utf-8", errors="replace") as fh:
                pka_text = fh.read()

        shifts = _parse_propka_output(pka_text)
        return True, shifts

    except FileNotFoundError:
        # propka3 not installed
        return False, []
    except asyncio.TimeoutError:
        return False, []
    except Exception:
        return False, []


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_catalytic_mechanism(
    gene: str,
    metal_ions: list[str],
    coordination_sites: list[MetalCoordination],
    pka_shifts: list[PKaShift],
) -> tuple[str, list[str], str]:
    """Call Gemini to synthesise catalytic mechanism + inhibition strategies.

    Returns (catalytic_mechanism, inhibition_strategies, synthesis_summary).
    """
    try:
        from core.gemini_interpreter import _call

        coord_residues = []
        for site in coordination_sites:
            coord_residues.append(
                f"{site.metal_ion} ({site.geometry}): {', '.join(site.coordinating_residues)}"
            )

        pka_shifts_data = [
            {
                "residue": s.residue,
                "model_pka": s.model_compound_pka,
                "predicted_pka": s.predicted_pka,
                "shift": s.shift,
                "significance": s.functional_significance,
            }
            for s in pka_shifts
        ]

        prompt = (
            f"Metalloenzyme: {gene}\n"
            f"Metal ions: {', '.join(metal_ions)}\n"
            f"Coordinating residues: {'; '.join(coord_residues) if coord_residues else 'Not determined'}\n"
            f"pKa-shifted residues: {_json.dumps(pka_shifts_data)}\n\n"
            "As a metalloenzyme expert, describe:\n"
            "1. The proposed catalytic mechanism\n"
            "2. How the metal coordination enables catalysis\n"
            "3. Three specific inhibition strategies (chelators, competitive inhibitors, allosteric)\n\n"
            "Return ONLY raw JSON (no markdown fences):\n"
            "{\n"
            '  "catalytic_mechanism": "...",\n'
            '  "inhibition_strategies": ["...", "...", "..."],\n'
            '  "synthesis": "2-3 sentence summary"\n'
            "}"
        )

        raw = await _call(prompt)
        if not raw:
            return "", [], ""

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        mechanism = str(data.get("catalytic_mechanism", "")).strip()
        strategies_raw = data.get("inhibition_strategies", [])
        strategies = [str(s).strip() for s in strategies_raw if s] if isinstance(strategies_raw, list) else []
        synthesis = str(data.get("synthesis", "")).strip()

        return mechanism, strategies, synthesis

    except Exception:
        return "", [], ""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_metalloenzyme_analysis(
    gene: str,
    cofactors: list[str],
    uniprot_entry: dict,
    af_pdb_url: Optional[str] = None,
    step_cb=None,
) -> MetalloenzymeAnalysis:
    """Analyse metal cofactor coordination, pKa perturbations, and catalytic mechanism.

    Auto-triggers when the cofactors list contains metal ions.

    Args:
        gene:          Gene symbol (e.g. "CA2", "SOD1").
        cofactors:     List of cofactor strings from UniProt extraction.
        uniprot_entry: Raw UniProt entry dict (used for binding-site features).
        af_pdb_url:    AlphaFold PDB download URL (optional; used for PROPKA).
        step_cb:       Optional async progress callback(msg: str).
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    result = MetalloenzymeAnalysis(gene=gene.upper())

    # ── 1. Detect metal cofactors ─────────────────────────────────────────────
    await _step("[1/6] Detecting metal cofactors…")
    metal_ions = _detect_metals(cofactors)
    if not metal_ions:
        # No metals detected — return empty result gracefully
        return result
    result.metal_ions = metal_ions

    # ── 2. Parse metal coordination residues from UniProt features ────────────
    await _step("[2/6] Parsing metal coordination residues from UniProt features…")
    try:
        coordination_sites = _parse_uniprot_binding_sites(uniprot_entry, metal_ions)
        result.coordination_sites = coordination_sites
    except Exception:
        coordination_sites = []

    # ── 3. Download AlphaFold structure for pKa analysis ─────────────────────
    await _step("[3/6] Downloading AlphaFold structure for pKa analysis…")
    pdb_tmp_path: Optional[str] = None
    if af_pdb_url:
        try:
            import httpx
            async with httpx.AsyncClient(timeout=60.0) as client:
                resp = await client.get(af_pdb_url)
                resp.raise_for_status()
                # Write to a named temp file that persists until we clean up
                with tempfile.NamedTemporaryFile(
                    suffix=".pdb", delete=False, mode="wb"
                ) as tmp:
                    tmp.write(resp.content)
                    pdb_tmp_path = tmp.name
        except Exception:
            pdb_tmp_path = None

    # ── 4. Run PROPKA pKa prediction ─────────────────────────────────────────
    await _step("[4/6] Running PROPKA pKa prediction…")
    pka_shifts: list[PKaShift] = []
    propka_available = False
    if pdb_tmp_path:
        try:
            propka_available, pka_shifts = await _run_propka(pdb_tmp_path)
        except Exception:
            propka_available = False
            pka_shifts = []
        finally:
            # Clean up temp file
            try:
                import os
                os.unlink(pdb_tmp_path)
            except Exception:
                pass
    result.propka_available = propka_available
    result.pka_shifts = pka_shifts

    # ── 5. Gemini catalytic mechanism synthesis ───────────────────────────────
    await _step("[5/6] Running Gemini catalytic mechanism synthesis…")
    try:
        mechanism, strategies, synthesis = await _synthesize_catalytic_mechanism(
            gene=gene,
            metal_ions=metal_ions,
            coordination_sites=coordination_sites,
            pka_shifts=pka_shifts,
        )
        result.catalytic_mechanism = mechanism
        result.inhibition_strategies = strategies
        result.gemini_synthesis = synthesis
    except Exception:
        pass

    # ── 6. Assemble metalloenzyme analysis ────────────────────────────────────
    await _step("[6/6] Assembling metalloenzyme analysis…")
    return result
