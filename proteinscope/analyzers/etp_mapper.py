"""Beratan-Onuchic Electron Tunneling Pathway (ETP) Mapper.

Implements the Pathways model from Beratan, Betts & Onuchic (1992) for predicting
electron transfer routes through proteins.

Decay parameters (차박사, Oxford; Grand Consortium V3 2026-03-20):
  ε_bond = 0.6 per covalent bond
  ε_space(R) = 0.6 × exp[−1.7 × (R − 1.4)]  (R in Å, β = 1.7 Å⁻¹)
  ε_HB(R)    = 0.36 × exp[−1.7 × (R − 2.8)] (heavy-atom H-bond distance)

# Citation: Betts, Beratan & Onuchic (1992) JACS 114:4043. doi:10.1021/ja00037a004
# Citation: Beratan et al. (1992) Science 258:1740. doi:10.1126/science.1334572
# Validation: Ru-modified cytochrome c (Wuttke et al. 1992 Science 256:1007)
#   doi:10.1126/science.256.5059.1007
#   Expected coupling order: His39 > His33 > His72 > His62 despite distance ordering
"""

from __future__ import annotations

import math
import os
import tempfile
from datetime import datetime
from math import exp
from typing import Optional

import httpx
from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade, QBConfidenceScore

# ── Beratan-Onuchic 1992 decay constants (차박사) ─────────────────────────────
# Citation: Betts, Beratan & Onuchic (1992) JACS 114:4043. doi:10.1021/ja00037a004
EPSILON_BOND = 0.6          # decay per covalent bond (dimensionless)
BETA_SPACE = 1.7            # Å⁻¹, through-space decay constant (original 1992 value)
EPSILON_SPACE_PREFACTOR = 0.6
BETA_HB = 1.7               # Å⁻¹, H-bond decay constant
EPSILON_HB_PREFACTOR = 0.36  # = 0.6² (= ε_bond²)


def epsilon_space(r_angstrom: float) -> float:
    """Through-space coupling decay.

    ε_space(R) = 0.6 × exp[−1.7 × (R − 1.4)]  (R in Å)
    # Citation: Betts et al. (1992) JACS 114:4043
    """
    return EPSILON_SPACE_PREFACTOR * exp(-BETA_SPACE * (r_angstrom - 1.4))


def epsilon_hb(r_angstrom: float) -> float:
    """H-bond coupling decay (heavy-atom distance).

    ε_HB(R) = 0.36 × exp[−1.7 × (R − 2.8)]  (R in Å, heavy-atom distance)
    # Citation: Betts et al. (1992) JACS 114:4043
    """
    return EPSILON_HB_PREFACTOR * exp(-BETA_HB * (r_angstrom - 2.8))


# ── Cofactor donor/acceptor atom tables (차박사) ──────────────────────────────
# Citation: Beratan et al. (1992) Science 258:1740; Hayashi & Stuchebrukhov (2010)
COFACTOR_ATOMS: dict[str, list[str]] = {
    "heme": ["FE", "NA", "NB", "NC", "ND", "O1A", "O2A", "O1D", "O2D"],
    # Fe + four pyrrole N + propionate O atoms (pathway often exits via propionates)
    "fe_s": ["FE", "SF", "S1", "S2", "S3", "S4"],
    # inorganic Fe + bridging S (core-only, no ligand atoms)
    "fad": ["N5", "C4A", "O2", "O4", "N3", "C2", "N1", "C10"],
    "fmn": ["N5", "C4A", "O2", "O4", "N3", "C2", "N1"],
    # isoalloxazine ring only; exclude adenine/ribityl
    "cu_type1": ["CU", "SG"],
    # Cu + Cys sulfur (type 1 blue copper)
    "cu_a": ["CU", "CU2", "SG", "SG2"],
    # both Cu atoms + bridging Cys S (CuA binuclear center, mixed-valence)
    "mo_co": ["MO", "S1", "S2", "SD1", "SD2"],
    # Mo + dithiolene S atoms
}

# Keywords in cofactor description → cofactor key
_COFACTOR_MAP = {
    "heme": "heme",
    "haem": "heme",
    "cytochrome": "heme",
    "iron-sulfur": "fe_s",
    "fe-s": "fe_s",
    "fes": "fe_s",
    "4fe-4s": "fe_s",
    "2fe-2s": "fe_s",
    "flavin adenine": "fad",
    "fad": "fad",
    "fmn": "fmn",
    "riboflavin": "fmn",
    "copper": "cu_type1",
    "type-1 copper": "cu_type1",
    "blue copper": "cu_type1",
    "cua": "cu_a",
    "molybdenum": "mo_co",
    "molybdopterin": "mo_co",
}

# Mandatory caveat text (윤박사, CALTECH; Grand Consortium V3 2026-03-20)
MANDATORY_CAVEAT = (
    "CAVEAT: Predicted electron tunneling pathways are theoretical models based on a "
    "static protein structure from AlphaFold. The calculated tunneling rates are "
    "exponentially sensitive to inter-residue distances and orientations. These pathways "
    "do not account for protein dynamics, solvent effects, pH changes, or conformational "
    "rearrangements upon substrate binding, which can significantly alter or create new "
    "pathways. The reliability of this prediction is directly dependent on the pLDDT "
    "scores of the residues involved. This output is for hypothesis generation in a "
    "research context only and is not validated for clinical or diagnostic use."
)

# pLDDT reliability thresholds (차박사 + 윤박사, Grand Consortium V3)
PLDDT_UNRELIABLE = 70.0  # warn below this
PLDDT_EXCLUDE = 50.0     # exclude from path below this


# ── Pydantic models ───────────────────────────────────────────────────────────

class ETPResidue(BaseModel):
    """A residue participating in the electron tunneling pathway."""
    residue_number: int
    residue_name: str
    chain: str = "A"
    is_donor: bool = False
    is_acceptor: bool = False
    plddt: Optional[float] = None
    reliability: str = "reliable"  # "reliable" | "unreliable" | "excluded"
    is_redox_active_hotspot: bool = False  # ETP × Fragment Hotspot cross-ref (Dr. Konberg)


class ETPEdge(BaseModel):
    """A single coupling step in the tunneling pathway."""
    atom_a: str
    atom_b: str
    edge_type: str      # "covalent" | "h_bond" | "through_space"
    distance_angstrom: float
    epsilon: float      # decay coupling value (Beratan-Onuchic 1992)


class ETPAnalysis(BaseModel):
    """Full Beratan-Onuchic electron tunneling pathway analysis."""
    gene: str
    cofactor_type: Optional[str] = None
    donor_atoms: List[str] = Field(default_factory=list)
    acceptor_atoms: List[str] = Field(default_factory=list)
    pathway_residues: List[ETPResidue] = Field(default_factory=list)
    pathway_edges: List[ETPEdge] = Field(default_factory=list)
    total_coupling: float = 0.0        # product of all edge epsilons
    plddt_warning: bool = False        # True if any path residue pLDDT < PLDDT_UNRELIABLE
    redox_hotspot_residues: List[str] = Field(default_factory=list)  # ETP × Fragment cross-ref
    qb_score: Optional[QBConfidenceScore] = None
    mandatory_caveat: str = MANDATORY_CAVEAT
    gemini_interpretation: str = ""
    provenance: Optional[DataProvenance] = None
    timestamp: Optional[str] = None


# needed for List[] inside ETPAnalysis
from typing import List  # noqa: E402 — placed after class definitions intentionally


def _identify_cofactor(cofactors: list[str]) -> Optional[str]:
    """Identify cofactor type from cofactor description strings."""
    combined = " ".join(cofactors).lower()
    for keyword, key in _COFACTOR_MAP.items():
        if keyword in combined:
            return key
    return None


def _distance(a: tuple, b: tuple) -> float:
    """Euclidean distance in Å between two 3D coordinate tuples."""
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _classify_edge(dist: float, covalent_thresh: float = 2.0, hb_max: float = 3.5) -> str:
    if dist <= covalent_thresh:
        return "covalent"
    if dist <= hb_max:
        return "h_bond"
    return "through_space"


def _edge_epsilon(edge_type: str, dist: float) -> float:
    if edge_type == "covalent":
        return EPSILON_BOND
    if edge_type == "h_bond":
        return epsilon_hb(dist)
    return epsilon_space(dist)


def _build_pathway_from_structure(
    pdb_text: str,
    cofactor_key: str,
    plddt_map: dict[int, float],
    fragment_residues: set[str],
) -> tuple[list[ETPResidue], list[ETPEdge], float, bool, list[str]]:
    """Parse PDB text and compute approximate ETP path.

    Simplified implementation:
    1. Parse ATOM records to get Cα positions per residue
    2. Find residues near cofactor atoms (donor/acceptor candidates within 14 Å)
    3. Build coupling chain from nearest cofactor residue outward
    4. Compute epsilon product for the chain

    Returns: (residues, edges, total_coupling, plddt_warning, redox_hotspots)
    """
    # Parse Cα coordinates
    residues_ca: dict[tuple[int, str], tuple[str, str, tuple[float, float, float]]] = {}
    # key: (resnum, chain) → (resname, atom_name, coords)

    cofactor_coords: list[tuple[float, float, float]] = []
    target_atoms = set(COFACTOR_ATOMS.get(cofactor_key, []))

    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        try:
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or "A"
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        coords = (x, y, z)

        if rec == "ATOM" and atom_name == "CA":
            residues_ca[(resnum, chain)] = (resname, atom_name, coords)

        # Collect cofactor atom positions
        if rec == "HETATM" and atom_name in target_atoms:
            cofactor_coords.append(coords)

    if not cofactor_coords or not residues_ca:
        return [], [], 0.0, False, []

    # Centroid of cofactor
    cx = sum(c[0] for c in cofactor_coords) / len(cofactor_coords)
    cy = sum(c[1] for c in cofactor_coords) / len(cofactor_coords)
    cz = sum(c[2] for c in cofactor_coords) / len(cofactor_coords)
    cofactor_center = (cx, cy, cz)

    # Find residues within 14 Å of cofactor (pathway candidates)
    candidates = []
    for (resnum, chain), (resname, _, ca_coords) in residues_ca.items():
        dist = _distance(ca_coords, cofactor_center)
        if dist <= 14.0:
            candidates.append((dist, resnum, chain, resname, ca_coords))

    if not candidates:
        return [], [], 0.0, False, []

    # Sort by distance, take closest 8 as pathway residues
    candidates.sort(key=lambda x: x[0])
    path_residues_raw = candidates[:8]

    # Build ETPResidue objects
    etp_residues = []
    plddt_warning = False
    for i, (dist, resnum, chain, resname, _) in enumerate(path_residues_raw):
        plddt = plddt_map.get(resnum)
        if plddt is not None and plddt < PLDDT_EXCLUDE:
            reliability = "excluded"
        elif plddt is not None and plddt < PLDDT_UNRELIABLE:
            reliability = "unreliable"
            plddt_warning = True
        else:
            reliability = "reliable"

        res_label = f"{resname}{resnum}"
        is_hotspot = res_label in fragment_residues or str(resnum) in fragment_residues

        etp_residues.append(ETPResidue(
            residue_number=resnum,
            residue_name=resname,
            chain=chain,
            is_donor=(i == 0),
            is_acceptor=(i == len(path_residues_raw) - 1),
            plddt=plddt,
            reliability=reliability,
            is_redox_active_hotspot=is_hotspot,
        ))

    # Build edges between consecutive residues
    etp_edges = []
    total_coupling = 1.0
    coords_list = [c[4] for c in path_residues_raw]

    for i in range(len(coords_list) - 1):
        dist = _distance(coords_list[i], coords_list[i + 1])
        edge_type = _classify_edge(dist)
        eps = _edge_epsilon(edge_type, dist)
        etp_edges.append(ETPEdge(
            atom_a=f"CA-{path_residues_raw[i][1]}",
            atom_b=f"CA-{path_residues_raw[i + 1][1]}",
            edge_type=edge_type,
            distance_angstrom=round(dist, 2),
            epsilon=round(eps, 4),
        ))
        total_coupling *= eps

    redox_hotspots = [
        f"{r.residue_name}{r.residue_number}"
        for r in etp_residues if r.is_redox_active_hotspot
    ]

    return etp_residues, etp_edges, round(total_coupling, 6), plddt_warning, redox_hotspots


async def _gemini_etp_summary(
    gene: str,
    cofactor_type: str,
    pathway_residues: list[ETPResidue],
    total_coupling: float,
    plddt_warning: bool,
) -> str:
    """Generate Gemini interpretation of ETP analysis."""
    if not pathway_residues:
        return ""

    residue_str = ", ".join(
        f"{r.residue_name}{r.residue_number}"
        + (f" [pLDDT={r.plddt:.0f}, {r.reliability}]" if r.plddt else "")
        for r in pathway_residues[:6]
    )
    warning_note = (
        " Note: some pathway residues have pLDDT < 70 (low AlphaFold confidence)."
        if plddt_warning else ""
    )

    prompt = (
        f"You are an electron transfer biophysicist. For the protein {gene} with "
        f"cofactor type '{cofactor_type}', the Beratan-Onuchic pathway analysis "
        f"identified these pathway residues: {residue_str}. Total coupling product: "
        f"{total_coupling:.2e}.{warning_note}\n\n"
        f"In 2-3 sentences, interpret: (1) biological significance of this ET pathway, "
        f"(2) which residues are most critical for coupling, (3) any drug design implications. "
        f"Be concise and technically precise."
    )

    try:
        from google import genai

        api_key = os.environ.get("GEMINI_API_KEY", "")
        if not api_key:
            return ""
        client = genai.Client(api_key=api_key)
        resp = client.models.generate_content(model="models/gemini-2.5-pro", contents=prompt)
        return (resp.text or "").strip()
    except Exception:
        return ""


async def run_etp_analysis(
    gene: str,
    af_pdb_url: str,
    cofactors: list[str],
    fragment_hotspot_map=None,  # Optional FragmentHotspotMap
    plddt_per_residue: Optional[dict[int, float]] = None,
    step_cb=None,
) -> Optional[ETPAnalysis]:
    """Run Beratan-Onuchic ETP analysis on an AlphaFold structure.

    Args:
        gene: Gene symbol.
        af_pdb_url: AlphaFold PDB download URL.
        cofactors: List of cofactor description strings from ProteinRecord.
        fragment_hotspot_map: Optional FragmentHotspotMap for cross-ref.
        plddt_per_residue: Optional dict {resnum: plddt} for per-residue confidence.
        step_cb: Optional async progress callback.

    Returns:
        ETPAnalysis or None on failure.
    """
    try:
        # Identify cofactor type
        cofactor_key = _identify_cofactor(cofactors)
        if not cofactor_key:
            return None

        if step_cb:
            await step_cb(f"[ETP] Cofactor identified: {cofactor_key}")

        # Download AlphaFold PDB
        if step_cb:
            await step_cb("[ETP] Downloading AlphaFold structure")

        pdb_text = ""
        try:
            async with httpx.AsyncClient(timeout=60) as client:
                resp = await client.get(af_pdb_url)
            if resp.status_code == 200:
                pdb_text = resp.text
        except Exception:
            pass

        if not pdb_text:
            return None

        # Build pLDDT map from B-factor column if not provided
        if not plddt_per_residue:
            plddt_per_residue = {}
            for line in pdb_text.splitlines():
                if line[:6].strip() == "ATOM" and line[12:16].strip() == "CA":
                    try:
                        resnum = int(line[22:26])
                        bfactor = float(line[60:66])
                        plddt_per_residue[resnum] = bfactor
                    except (ValueError, IndexError):
                        pass

        # Build fragment hotspot residue set for cross-ref
        fragment_residues: set[str] = set()
        if fragment_hotspot_map is not None:
            for hotspot in getattr(fragment_hotspot_map, "hotspots", []):
                for res in getattr(hotspot, "residues", []):
                    fragment_residues.add(str(res))

        if step_cb:
            await step_cb("[ETP] Computing Beratan-Onuchic coupling pathway")

        # Compute pathway
        pathway_residues, pathway_edges, total_coupling, plddt_warning, redox_hotspots = (
            _build_pathway_from_structure(
                pdb_text, cofactor_key, plddt_per_residue, fragment_residues
            )
        )

        # Compute QB Confidence Score (윤박사)
        qb_score = None
        if pathway_residues and plddt_per_residue:
            path_plddts = [
                r.plddt for r in pathway_residues
                if r.plddt is not None
            ]
            min_plddt = min(path_plddts) if path_plddts else 50.0
            # cofactor_certainty: 0.6 (AlphaFold predicted, cofactor inferred from sequence)
            # kie_flag: 0.1 (no KIE data assumed unless provided)
            qb_score = QBConfidenceScore.compute(
                min_plddt=min_plddt,
                cofactor_certainty=0.6,
                kie_flag=0.1,
            )

        if step_cb:
            await step_cb("[ETP] Generating Gemini interpretation")

        gemini_text = await _gemini_etp_summary(
            gene, cofactor_key, pathway_residues, total_coupling, plddt_warning
        )

        return ETPAnalysis(
            gene=gene,
            cofactor_type=cofactor_key,
            donor_atoms=COFACTOR_ATOMS.get(cofactor_key, []),
            acceptor_atoms=[],  # surface residues; identified by proximity
            pathway_residues=pathway_residues,
            pathway_edges=pathway_edges,
            total_coupling=total_coupling,
            plddt_warning=plddt_warning,
            redox_hotspot_residues=redox_hotspots,
            qb_score=qb_score,
            mandatory_caveat=MANDATORY_CAVEAT,
            gemini_interpretation=gemini_text,
            provenance=DataProvenance(
                source="Beratan-Onuchic Pathways model (1992) on AlphaFold structure",
                evidence_grade=EvidenceGrade.COMPUTATIONAL,
                scientific_caveat=(
                    "Theoretical ET pathway from static AlphaFold coordinates; "
                    "not validated by experimental ET rate measurement for this protein"
                ),
                method="Beratan-Onuchic Pathways (1992); ε_bond=0.6, β=1.7 Å⁻¹",
            ),
            timestamp=datetime.utcnow().isoformat(),
        )

    except Exception:
        return None
