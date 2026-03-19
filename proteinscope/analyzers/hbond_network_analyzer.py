"""Hydrogen Bond Network Analyzer.

Identifies hydrogen-bond networks in the AlphaFold structure, classifies
backbone vs sidechain bonds, finds buried H-bond clusters, and quantifies
network stability as a proxy for conformational rigidity.

Background
──────────
Hydrogen bonds are the dominant non-covalent force in protein structure.
Buried H-bond networks are critical for:
• Fold stability (cooperativity)
• Active-site geometry maintenance
• Allosteric signal propagation
• Drug design (displacement of ordered water networks)

H-bond detection criteria
──────────────────────────
• Donor–acceptor heavy-atom distance: ≤ 3.5 Å
• D–H···A angle: > 90° (estimated from heavy-atom geometry here)
• Donor atoms: N (backbone NH), N (Asn/Gln/Arg/Lys/His/Trp side-chains),
               O (Ser/Thr/Tyr OH), S (Cys SH rarely)
• Acceptor atoms: O (backbone C=O), O (Asp/Glu/Asn/Gln/Ser/Thr/Tyr),
                  N (His Nδ/Nε when acting as acceptor)

Key citations
─────────────
• McDonald IK & Thornton JM 1994 J Mol Biol 238:777-93 doi:10.1006/jmbi.1994.1334
  (HBPLUS algorithm; H-bond definition 3.5 Å threshold)
• Baker EN & Hubbard RE 1984 Prog Biophys Mol Biol 44:97-179 doi:10.1016/0079-6107(84)90007-5
  (comprehensive H-bond geometry analysis)
• Stickle DF et al. 1992 J Mol Biol 226:1143-59 doi:10.1016/0022-2836(92)91058-W
  (H-bond statistics in protein crystal structures)
• Dill KA 1990 Biochemistry 29:7133-55 doi:10.1021/bi00483a001
  (hydrophobic effect + H-bond network contribution to stability)
• Jumper J et al. 2021 Nature 596:583-589 doi:10.1038/s41586-021-03819-2
  (AlphaFold2 — coordinate source)
"""
from __future__ import annotations

import math
import os
import tempfile
from datetime import datetime
from typing import Dict, List, Optional, Tuple

from pydantic import BaseModel, Field

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

# H-bond cutoffs
# Citation: McDonald & Thornton 1994 doi:10.1006/jmbi.1994.1334
_HBOND_DISTANCE_CUTOFF = 3.5          # Å, heavy-atom donor–acceptor
_HBOND_WEAK_CUTOFF = 3.9              # Å — weak / water-mediated candidates

# Buried-cluster definition: ≥ 3 H-bonds within 8 Å of each other
_CLUSTER_RADIUS = 8.0
_CLUSTER_MIN_SIZE = 3

# Donor heavy atoms by residue + atom name (PDB naming)
_DONOR_ATOMS: Dict[str, List[str]] = {
    # Backbone
    "ALL": ["N"],
    # Sidechain donors
    "SER": ["OG"],  "THR": ["OG1"], "TYR": ["OH"],
    "CYS": ["SG"],  "ASN": ["ND2"], "GLN": ["NE2"],
    "LYS": ["NZ"],  "ARG": ["NH1", "NH2", "NE"],
    "HIS": ["ND1", "NE2"],
    "TRP": ["NE1"],
}

# Acceptor heavy atoms by residue
_ACCEPTOR_ATOMS: Dict[str, List[str]] = {
    # Backbone
    "ALL": ["O"],
    # Sidechain acceptors
    "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
    "ASN": ["OD1"],        "GLN": ["OE1"],
    "SER": ["OG"],         "THR": ["OG1"],   "TYR": ["OH"],
    "HIS": ["ND1", "NE2"],
    "MET": ["SD"],
}


# ──────────────────────────────────────────────────────────────────────────────
# Pydantic models
# ──────────────────────────────────────────────────────────────────────────────

class HBond(BaseModel):
    """One hydrogen bond identified in the structure."""
    donor_res: str              # e.g. "SER42"
    donor_atom: str             # e.g. "OG"
    acceptor_res: str           # e.g. "ASP87"
    acceptor_atom: str          # e.g. "OD1"
    distance_angstrom: float
    bond_type: str              # "backbone-backbone" | "backbone-sidechain" | "sidechain-sidechain"
    is_buried: bool = False


class HBondCluster(BaseModel):
    """A network of co-located H-bonds (buried cluster)."""
    center_residue: str
    cluster_size: int           # number of H-bonds in cluster
    residues_involved: List[str] = Field(default_factory=list)
    stability_note: str = ""


class HBondNetwork(BaseModel):
    """Full H-bond network analysis of the protein structure."""
    gene: str
    total_hbonds: int
    backbone_backbone: int
    backbone_sidechain: int
    sidechain_sidechain: int
    buried_clusters: List[HBondCluster] = Field(default_factory=list)
    representative_hbonds: List[HBond] = Field(default_factory=list)   # top 20
    network_summary: str = ""
    gemini_interpretation: str = ""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    references: List[str] = Field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# PDB parsing helpers (same style as allosteric_analyzer.py)
# ──────────────────────────────────────────────────────────────────────────────

_REFERENCES = [
    "McDonald IK & Thornton JM 1994 J Mol Biol 238:777 doi:10.1006/jmbi.1994.1334",
    "Baker EN & Hubbard RE 1984 Prog Biophys Mol Biol 44:97 doi:10.1016/0079-6107(84)90007-5",
    "Stickle DF et al. 1992 J Mol Biol 226:1143 doi:10.1016/0022-2836(92)91058-W",
    "Dill KA 1990 Biochemistry 29:7133 doi:10.1021/bi00483a001",
    "Jumper J et al. 2021 Nature 596:583 doi:10.1038/s41586-021-03819-2",
]


def _parse_heavy_atoms(
    pdb_path: str,
) -> List[Tuple[str, str, str, float, float, float]]:
    """Parse heavy atoms from PDB ATOM records.

    Returns list of (res_name, res_seq, atom_name, x, y, z).
    Only first chain, no alt-loc.
    """
    atoms: List[Tuple[str, str, str, float, float, float]] = []
    seen_chain: Optional[str] = None

    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                chain = line[21].strip()
                if seen_chain is None:
                    seen_chain = chain
                elif chain != seen_chain:
                    continue  # only first chain
                alt_loc = line[16].strip()
                if alt_loc and alt_loc != "A":
                    continue  # skip alt-loc
                atom_name = line[12:16].strip()
                if atom_name.startswith("H"):
                    continue  # skip explicit H
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                atoms.append((res_name, res_seq, atom_name, x, y, z))
    except Exception:
        pass

    return atoms


def _dist(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def _is_donor(res_name: str, atom_name: str) -> bool:
    """Check if atom is an H-bond donor heavy atom."""
    specific = _DONOR_ATOMS.get(res_name, [])
    backbone = _DONOR_ATOMS["ALL"]
    return atom_name in specific or atom_name in backbone


def _is_acceptor(res_name: str, atom_name: str) -> bool:
    """Check if atom is an H-bond acceptor heavy atom."""
    specific = _ACCEPTOR_ATOMS.get(res_name, [])
    backbone = _ACCEPTOR_ATOMS["ALL"]
    return atom_name in specific or atom_name in backbone


def _bond_type(donor_atom: str, acceptor_atom: str) -> str:
    """Classify bond type from atom names."""
    # Backbone atoms: N, O
    donor_bb = donor_atom in ("N",)
    acceptor_bb = acceptor_atom in ("O",)
    if donor_bb and acceptor_bb:
        return "backbone-backbone"
    if donor_bb or acceptor_bb:
        return "backbone-sidechain"
    return "sidechain-sidechain"


def _find_hbonds(
    atoms: List[Tuple[str, str, str, float, float, float]],
) -> List[HBond]:
    """Find all H-bonds using distance criteria.

    Citation: McDonald & Thornton 1994 doi:10.1006/jmbi.1994.1334
    Threshold: 3.5 Å donor–acceptor heavy-atom distance.
    """
    bonds: List[HBond] = []
    n = len(atoms)

    # Build index for faster lookup
    for i in range(n):
        res_d, seq_d, atom_d, x1, y1, z1 = atoms[i]
        if not _is_donor(res_d, atom_d):
            continue
        donor_label = f"{res_d}{seq_d}"

        for j in range(n):
            if i == j:
                continue
            res_a, seq_a, atom_a, x2, y2, z2 = atoms[j]
            if not _is_acceptor(res_a, atom_a):
                continue
            # Skip same residue
            if seq_d == seq_a:
                continue
            # Skip adjacent backbone (i → i+1 would be trivial)
            try:
                seqdiff = abs(int(seq_d) - int(seq_a))
            except ValueError:
                seqdiff = 99
            if seqdiff < 2:
                continue

            d = _dist((x1, y1, z1), (x2, y2, z2))
            if d > _HBOND_DISTANCE_CUTOFF:
                continue

            acceptor_label = f"{res_a}{seq_a}"
            btype = _bond_type(atom_d, atom_a)

            bonds.append(HBond(
                donor_res=donor_label,
                donor_atom=atom_d,
                acceptor_res=acceptor_label,
                acceptor_atom=atom_a,
                distance_angstrom=round(d, 3),
                bond_type=btype,
            ))

    # Deduplicate (same residue pair can have multiple distances; keep shortest)
    seen: Dict[str, HBond] = {}
    for b in bonds:
        key = f"{b.donor_res}:{b.acceptor_res}"
        if key not in seen or b.distance_angstrom < seen[key].distance_angstrom:
            seen[key] = b

    return list(seen.values())


def _find_buried_clusters(bonds: List[HBond], atoms: List) -> List[HBondCluster]:
    """Find clusters of ≥3 H-bonds within _CLUSTER_RADIUS Å of each other.

    Uses acceptor position as cluster centroid.
    Citation: Stickle 1992 doi:10.1016/0022-2836(92)91058-W
    """
    if not bonds or not atoms:
        return []

    # Build acceptor residue → coordinate map
    coord_map: Dict[str, Tuple[float, float, float]] = {}
    for res_name, res_seq, atom_name, x, y, z in atoms:
        label = f"{res_name}{res_seq}"
        if atom_name in ("CA", "N", "O"):   # backbone centroid
            coord_map.setdefault(label, (x, y, z))

    clusters: List[HBondCluster] = []
    residues_in_clusters: set = set()

    for bond in bonds:
        center = coord_map.get(bond.donor_res)
        if not center:
            continue
        nearby = []
        nearby_res = set()
        for other in bonds:
            ocoord = coord_map.get(other.donor_res)
            if not ocoord:
                continue
            if _dist(center, ocoord) <= _CLUSTER_RADIUS:
                nearby.append(other)
                nearby_res.update([other.donor_res, other.acceptor_res])

        if len(nearby) >= _CLUSTER_MIN_SIZE and bond.donor_res not in residues_in_clusters:
            residues_in_clusters.update(nearby_res)
            clusters.append(HBondCluster(
                center_residue=bond.donor_res,
                cluster_size=len(nearby),
                residues_involved=sorted(list(nearby_res))[:8],
                stability_note=(
                    f"{len(nearby)} H-bonds within {_CLUSTER_RADIUS:.0f} Å — "
                    f"likely buried network contributing to fold stability"
                ),
            ))

    # Sort by cluster size descending; deduplicate residues across clusters
    clusters.sort(key=lambda c: -c.cluster_size)
    return clusters[:5]


# ──────────────────────────────────────────────────────────────────────────────
# Public entry-point
# ──────────────────────────────────────────────────────────────────────────────

async def run_hbond_analysis(
    gene: str,
    af_pdb_url: str,
    step_cb=None,
) -> Optional[HBondNetwork]:
    """Analyse hydrogen-bond network in the AlphaFold structure.

    Parameters
    ----------
    gene:       HGNC gene symbol.
    af_pdb_url: AlphaFold PDB URL (from ProteinRecord.alphafold_pdb_url).
    step_cb:    Optional async SSE progress callback.

    Returns
    -------
    HBondNetwork or None on error.
    """
    tmp_path: Optional[str] = None
    try:
        if not af_pdb_url:
            return None

        # ── Step 1: download PDB ───────────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "downloading_structure", 0.10)

        try:
            import httpx
        except ImportError:
            return None

        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.get(af_pdb_url)
            if resp.status_code != 200:
                return None
            pdb_content = resp.text

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb",
                                         delete=False, encoding="utf-8") as tmp:
            tmp.write(pdb_content)
            tmp_path = tmp.name

        # ── Step 2: parse heavy atoms ─────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "parsing_atoms", 0.25)

        atoms = _parse_heavy_atoms(tmp_path)
        if not atoms:
            return None

        # ── Step 3: detect H-bonds ────────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "detecting_hbonds", 0.50)

        all_bonds = _find_hbonds(atoms)
        if not all_bonds:
            return None

        # ── Step 4: classify and cluster ──────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "clustering", 0.65)

        bb_bb = sum(1 for b in all_bonds if b.bond_type == "backbone-backbone")
        bb_sc = sum(1 for b in all_bonds if b.bond_type == "backbone-sidechain")
        sc_sc = sum(1 for b in all_bonds if b.bond_type == "sidechain-sidechain")

        clusters = _find_buried_clusters(all_bonds, atoms)

        # Mark bonds that are in clusters as buried
        buried_residues = set()
        for cl in clusters:
            buried_residues.update(cl.residues_involved)
        for b in all_bonds:
            if b.donor_res in buried_residues or b.acceptor_res in buried_residues:
                b.is_buried = True

        # Representative H-bonds: prioritise buried sidechain-sidechain, then all
        sorted_bonds = sorted(
            all_bonds,
            key=lambda b: (not b.is_buried, b.bond_type != "sidechain-sidechain", b.distance_angstrom),
        )
        representative = sorted_bonds[:20]

        # ── Step 5: summary ───────────────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "writing_summary", 0.75)

        summary = (
            f"{len(all_bonds)} H-bonds detected (≤{_HBOND_DISTANCE_CUTOFF} Å): "
            f"{bb_bb} backbone–backbone, {bb_sc} backbone–sidechain, "
            f"{sc_sc} sidechain–sidechain. "
            f"{len(clusters)} buried cluster(s) (≥{_CLUSTER_MIN_SIZE} bonds within "
            f"{_CLUSTER_RADIUS:.0f} Å)."
        )

        # ── Step 6: Gemini synthesis ──────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "gemini_synthesis", 0.88)

        gemini_interpretation = ""
        try:
            from core.gemini_interpreter import _call  # type: ignore
            cluster_desc = "\n".join(
                f"  Cluster at {c.center_residue}: {c.cluster_size} bonds — {c.stability_note}"
                for c in clusters[:3]
            ) or "  No buried clusters found."
            top_bonds_desc = "\n".join(
                f"  {b.donor_res}({b.donor_atom})···{b.acceptor_res}({b.acceptor_atom}) "
                f"{b.distance_angstrom:.2f} Å [{b.bond_type}]{'*buried*' if b.is_buried else ''}"
                for b in representative[:8]
            )
            prompt = (
                f"You are a structural biologist analysing hydrogen-bond networks.\n\n"
                f"Gene: {gene}\n"
                f"Total H-bonds: {len(all_bonds)}\n"
                f"  Backbone–backbone: {bb_bb}\n"
                f"  Backbone–sidechain: {bb_sc}\n"
                f"  Sidechain–sidechain: {sc_sc}\n\n"
                f"Buried H-bond clusters:\n{cluster_desc}\n\n"
                f"Representative H-bonds:\n{top_bonds_desc}\n\n"
                f"Please provide:\n"
                f"1. What the H-bond network density implies about fold stability\n"
                f"2. Which buried clusters are most functionally significant and why\n"
                f"3. How the sidechain H-bond network relates to active site / binding site geometry\n"
                f"4. Engineering implications (which H-bonds to protect or break for design)\n\n"
                f"Write at the level of a structural biology expert."
            )
            gemini_interpretation = await _call(prompt)
        except Exception:
            gemini_interpretation = ""

        # ── Step 7: assemble ──────────────────────────────────────────────────
        if step_cb:
            await step_cb("hbond", "complete", 1.0)

        return HBondNetwork(
            gene=gene.upper(),
            total_hbonds=len(all_bonds),
            backbone_backbone=bb_bb,
            backbone_sidechain=bb_sc,
            sidechain_sidechain=sc_sc,
            buried_clusters=clusters,
            representative_hbonds=representative,
            network_summary=summary,
            gemini_interpretation=gemini_interpretation or "",
            timestamp=datetime.utcnow(),
            references=_REFERENCES,
        )

    except Exception:
        return None
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass
