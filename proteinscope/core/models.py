from __future__ import annotations

from pydantic import BaseModel, Field
from typing import Optional, List

from core.evidence import EvidenceGrade, DataProvenance  # noqa: F401 — re-exported


class SequenceFeature(BaseModel):
    feature_type: str               # e.g. "binding site", "active site", "domain"
    start: int
    end: int
    sequence_fragment: str          # extracted sub-sequence
    description: Optional[str] = None
    ligand: Optional[str] = None    # what molecule binds here


class Isoform(BaseModel):
    isoform_id: str
    name: str
    sequence: str
    differences: str                # description of how it differs from canonical


class CrossSpeciesEntry(BaseModel):
    organism: str
    uniprot_id: str
    identity_pct: float             # sequence identity to query protein
    interaction_domain_seq: Optional[str] = None
    compatible: Optional[bool] = None  # whether binding domains are conserved


class ClinicalEntry(BaseModel):
    variant: Optional[str] = None
    disease: str
    significance: str               # pathogenic / likely pathogenic / etc.
    omim_id: Optional[str] = None
    pubmed_ids: List[str] = Field(default_factory=list)
    provenance: Optional[DataProvenance] = None   # P0: evidence grading
    # P1: ESM-2 variant fitness fields (populated by query_engine after scoring)
    esm2_score: Optional[float] = None            # ΔLL; more negative = more damaging
    esm2_interpretation: Optional[str] = None     # "Benign" / "Possibly damaging" / "Likely damaging"


class Reference(BaseModel):
    pubmed_id: Optional[str] = None
    doi: Optional[str] = None
    title: str
    authors: List[str] = Field(default_factory=list)
    journal: str
    year: int
    apa_citation: str = ""
    vancouver_citation: str = ""
    bibtex: str = ""


class ProteinRecord(BaseModel):
    # --- Identity ---
    query: str
    uniprot_id: str
    protein_name: str
    gene_name: str
    organism: str
    taxonomy_id: str

    # --- Function ---
    function_description: str = ""
    biological_process: List[str] = Field(default_factory=list)   # GO terms
    molecular_function: List[str] = Field(default_factory=list)   # GO terms
    cofactors: List[str] = Field(default_factory=list)
    optimal_temperature_celsius: Optional[float] = None
    optimal_ph: Optional[float] = None
    subcellular_location: List[str] = Field(default_factory=list)

    # --- Sequences ---
    canonical_sequence: str = ""    # amino acid FASTA
    sequence_length: int = 0
    encoding_dna_cds: Optional[str] = None   # CDS nucleotide sequence
    is_from_spliced_mrna: bool = False
    isoforms: List[Isoform] = Field(default_factory=list)

    # --- Interaction Domains ---
    binding_sites: List[SequenceFeature] = Field(default_factory=list)
    active_sites: List[SequenceFeature] = Field(default_factory=list)
    domains: List[SequenceFeature] = Field(default_factory=list)

    # --- Clinical ---
    clinical_variants: List[ClinicalEntry] = Field(default_factory=list)
    deficiency_diseases: List[str] = Field(default_factory=list)

    # --- Structure ---
    alphafold_pdb_url: Optional[str] = None
    alphafold_plddt_score: Optional[float] = None   # overall confidence (0-100)
    experimental_pdb_ids: List[str] = Field(default_factory=list)

    # --- Cross-species ---
    cross_species_comparison: List[CrossSpeciesEntry] = Field(default_factory=list)

    # --- References ---
    references: List[Reference] = Field(default_factory=list)

    # --- AI summary (Gemini) ---
    ai_summary: Optional[str] = None

    # --- P1: ESM-2 variant fitness landscape ---
    variant_fitness_scores: List["VariantScore"] = Field(default_factory=list)

    # --- P2: PTM regulatory logic ---
    ptm_logic: Optional["PTMLogicGraph"] = None

    # --- Pathway extension fields ---
    drug_interactions: List["DrugInteraction"] = Field(default_factory=list)
    pgx_variants: List["PGxVariant"] = Field(default_factory=list)
    pharmacogenomic_pathways: List["PharmacogenomicPathway"] = Field(default_factory=list)
    metabolic_pathways: List["MetabolicPathway"] = Field(default_factory=list)
    signaling_pathways: List["SignalingPathway"] = Field(default_factory=list)
    protein_interactions: List["ProteinInteraction"] = Field(default_factory=list)
    disease_pathways: List["DiseasePathway"] = Field(default_factory=list)
    annotated_diagrams: List["AnnotatedDiagram"] = Field(default_factory=list)

    # --- P3: IDP/LLPS analysis (auto-triggers when pLDDT < 70) ---
    idp_analysis: Optional["IDPAnalysis"] = None

    # --- P4: Metalloenzyme / PROPKA (auto-triggers when metal cofactors detected) ---
    metalloenzyme_analysis: Optional["MetalloenzymeAnalysis"] = None

    # --- P5: Host-Protein Expression Compatibility Matrix ---
    expression_compatibility: Optional["ExpressionCompatibilityReport"] = None

    # --- P6: Binding Pocket Druggability ---
    druggability_report: Optional["DruggabilityReport"] = None

    # --- Ph2: Conformational ANM/GNM analysis ---
    conformational_analysis: Optional["ConformationalAnalysis"] = None

    # --- P9: Inverse Folding designs ---
    inverse_designs: List["InverseDesign"] = Field(default_factory=list)

    # --- P10: Therapeutic Decision Simulator ---
    therapeutic_decision: Optional["TherapeuticDecisionReport"] = None

    # --- P11: Antigen Discovery Mode ---
    antigen_discovery: Optional["AntigenDiscoveryReport"] = None

    # --- Layer 1: Molecular State ---
    allosteric_network: Optional["AllostericNetwork"] = None
    epistasis_report: Optional["EpistasisReport"] = None
    binding_energy_estimates: List["BindingEnergyEstimate"] = Field(default_factory=list)
    observables_prediction: Optional["ObservablesPrediction"] = None

    # --- Layer 2: Predictive Mechanism ---
    active_learning_recommendation: Optional["MutationRecommendation"] = None

    # --- Layer 3: Engineering Feasibility ---
    glycan_analysis: Optional["GlycanAnalysis"] = None
    admet_report: Optional["ADMETReport"] = None
    directed_evolution_plan: Optional["DirectedEvolutionPlan"] = None
    protac_feasibility: Optional["PRORTACFeasibilityReport"] = None
    selectivity_landscape: Optional["SelectivityLandscape"] = None
    fragment_hotspot_map: Optional["FragmentHotspotMap"] = None
    aav_design: Optional["AAVDesignReport"] = None
    crispr_plan: Optional["CRISPRIntegrationPlan"] = None
    genetic_circuit_design: Optional["GeneticCircuitDesign"] = None

    # --- Layer 3+: Additional Engineering ---
    covalent_inhibitor_design: Optional["CovalentInhibitorDesign"] = None
    engineering_strategy: Optional["EngineeringStrategyReport"] = None
    hbond_network: Optional["HBondNetwork"] = None

    # --- V3 Grand Consortium: Disease Association Layer ---
    disease_association_report: Optional["DiseaseAssociationReport"] = None

    # --- V3 Grand Consortium: ETP Mapper (Beratan-Onuchic) ---
    etp_analysis: Optional["ETPAnalysis"] = None


# ── Drug & Pharmacogenomics ──────────────────────────────────────────────────

class DrugInteraction(BaseModel):
    drug_name: str
    drug_id: str                    # ChEMBL ID or PharmGKB ID
    interaction_type: str           # inhibitor / activator / substrate / inducer
    mechanism: Optional[str] = None  # e.g. "competitive inhibition at active site"
    clinical_significance: Optional[str] = None  # PharmGKB evidence level 1A-4
    fda_label_mentioned: bool = False
    sources: List[str] = Field(default_factory=list)  # ["PharmGKB", "DGIdb", "ChEMBL"]
    provenance: Optional[DataProvenance] = None   # P0: evidence grading


class PGxVariant(BaseModel):
    variant_id: str                 # rsID or HGVS notation
    gene_symbol: str
    drug_name: str
    phenotype: str                  # e.g. "reduced metabolism", "toxicity risk"
    evidence_level: str             # PharmGKB 1A / 1B / 2A / 2B / 3 / 4
    population: Optional[str] = None  # ethnic group if specified
    pubmed_ids: List[str] = Field(default_factory=list)


class PharmacogenomicPathway(BaseModel):
    pathway_id: str                 # e.g. "PA150654557"
    pathway_name: str               # e.g. "Warfarin Pathway, Pharmacodynamics"
    drugs_involved: List[str] = Field(default_factory=list)
    genes_involved: List[str] = Field(default_factory=list)
    diagram_url: str
    diagram_image_path: Optional[str] = None  # local PNG for PDF embedding


# ── Metabolic ────────────────────────────────────────────────────────────────

class MetabolicReaction(BaseModel):
    reaction_id: str                # KEGG R-number
    equation: str                   # e.g. "ATP + Glucose -> ADP + Glucose-6P"
    substrates: List[str] = Field(default_factory=list)
    products: List[str] = Field(default_factory=list)
    cofactors_required: List[str] = Field(default_factory=list)
    enzyme_role: str = ""           # e.g. "catalyst", "allosteric activator"
    reversible: bool = False


class MetabolicPathway(BaseModel):
    pathway_id: str                 # KEGG map ID e.g. "hsa00010"
    pathway_name: str               # e.g. "Glycolysis / Gluconeogenesis"
    protein_role: str = ""          # role of the queried protein in this pathway
    reactions: List[MetabolicReaction] = Field(default_factory=list)
    upstream_proteins: List[str] = Field(default_factory=list)
    downstream_proteins: List[str] = Field(default_factory=list)
    diagram_url: str
    diagram_image_path: Optional[str] = None


# ── Signaling & Regulatory ───────────────────────────────────────────────────

class SignalingPathway(BaseModel):
    pathway_id: str                 # Reactome R-HSA-XXXXXXX
    pathway_name: str
    hierarchy: List[str] = Field(default_factory=list)  # breadcrumb
    protein_role: str = ""          # e.g. "kinase", "scaffold", "transcription factor"
    activates: List[str] = Field(default_factory=list)
    inhibits: List[str] = Field(default_factory=list)
    activated_by: List[str] = Field(default_factory=list)
    inhibited_by: List[str] = Field(default_factory=list)
    diagram_url: str
    diagram_image_path: Optional[str] = None


class ProteinInteraction(BaseModel):
    partner_protein: str
    partner_gene: str
    interaction_type: str           # "physical", "functional", "co-expression"
    confidence_score: float         # STRING score 0.0 - 1.0
    experimentally_confirmed: bool = False


# ── Disease Association Layer (V3 Grand Consortium) ──────────────────────────

class DiseaseAssociation(BaseModel):
    """Single gene-disease association from multi-DB integration.

    Sources: OpenTargets (GraphQL) + DisGeNET (REST) + ClinGen (REST) + HPO (REST)
    Dedup strategy: canonical_id uses EFO > MONDO > UMLS prefix priority.
    # Citation: 김박사 (Codex GPT-5.4) + 노박사 (Gemini 2.5 Pro), Grand Consortium V3 2026-03-20
    """
    disease_id: str                              # canonical ID (EFO > MONDO > UMLS)
    disease_name: str
    canonical_mondo_id: Optional[str] = None
    source_ids: List[str] = Field(default_factory=list)  # all IDs across sources
    score: float                                 # 0-1, OpenTargets overall score
    source: str                                  # "OpenTargets" | "DisGeNET" | "ClinGen"
    evidence_count: Optional[int] = None
    evidence_type: str = "genetic_association"   # genetic | somatic | drug_target | literature
    # LoF/GoF from OpenTargets variantFunctionalConsequenceId (노박사):
    # LoF: SO:0001587/0001589/0001574/0001575 | GoF: SO:0001583 + context
    lof_gof_hint: Optional[str] = None          # "LoF" | "GoF" | None
    therapeutic_modality_hint: Optional[str] = None  # "inhibitor" | "replacement" | "gene_therapy"
    hpo_terms: List[str] = Field(default_factory=list)
    clingen_classification: Optional[str] = None  # "Definitive" | "Strong" | "Moderate" | "Limited"
    mechanistic_summary: Optional[str] = None    # 노박사's template: Gene variant → mechanism → pathway → phenotype
    provenance: Optional[DataProvenance] = None


class DiseaseAssociationReport(BaseModel):
    gene: str
    associations: List[DiseaseAssociation] = Field(default_factory=list)
    total_count: int = 0
    gemini_summary: str = ""
    timestamp: Optional[str] = None


# ── Disease / Pathology ───────────────────────────────────────────────────────

class PathologyStep(BaseModel):
    step_number: int
    event: str                      # e.g. "Loss of CFTR function"
    consequence: str                # e.g. "Thick mucus accumulation in airways"
    reversible: bool = False
    therapeutic_target: bool = False  # is this step a known drug target?


class DiseasePathway(BaseModel):
    disease_name: str
    omim_id: Optional[str] = None
    pathway_id: str
    pathway_source: str             # "Reactome" or "KEGG"
    cascade: List[PathologyStep] = Field(default_factory=list)
    diagram_url: str
    diagram_image_path: Optional[str] = None


# ── Annotation Extension ──────────────────────────────────────────────────────

class NodeAnnotation(BaseModel):
    callout_number: int = 0
    node_id: str = ""
    gene_symbol: str
    protein_name: str
    x: int
    y: int
    node_width: int = 46
    node_height: int = 17
    arrow_label: str = ""
    interaction_comment: str = ""
    drug_comment: Optional[str] = None
    disease_comment: Optional[str] = None
    upstream_comment: Optional[str] = None
    downstream_comment: Optional[str] = None
    polished_comment: Optional[str] = None
    short_comment: str = ""
    extended_footnote: Optional[str] = None
    source: str = ""
    diagram_source_url: str = ""


class AnnotatedDiagram(BaseModel):
    pathway_id: str
    pathway_name: str
    source: str
    plain_image_path: str = ""
    annotated_image_path: str = ""
    external_url: str = ""
    annotations: List["NodeAnnotation"] = Field(default_factory=list)
    footnote_text: str = ""


# ── P1: ESM-2 Variant Fitness ─────────────────────────────────────────────────
# Re-export VariantScore from analyzers so callers only need to import from models.

try:
    from analyzers.variant_scorer import VariantScore  # noqa: F401
except Exception:
    # fair-esm not installed yet — define a minimal stub so ProteinRecord stays valid
    class VariantScore(BaseModel):  # type: ignore[no-redef]
        position: int
        ref_aa: str
        alt_aa: str
        variant: str
        delta_log_likelihood: Optional[float] = None
        wild_type_ll: Optional[float] = None
        alt_ll: Optional[float] = None
        predicted_impact: Optional[str] = None
        model_name: str = ""


# ── P2: PTM Logic Engine ──────────────────────────────────────────────────────
try:
    from analyzers.ptm_logic_engine import PTMSite, PTMLogicRelation, PTMLogicGraph  # noqa: F401
except Exception:
    pass

# ── P3: IDP / LLPS Analysis ───────────────────────────────────────────────────
try:
    from analyzers.idp_analyzer import IDPRegion, IDPAnalysis  # noqa: F401
except Exception:
    pass

# ── P4: Metalloenzyme Analysis ────────────────────────────────────────────────
try:
    from analyzers.metalloenzyme_analyzer import MetalCoordination, MetalloenzymeAnalysis  # noqa: F401
except Exception:
    pass

# ── P5: Expression Compatibility ──────────────────────────────────────────────
try:
    from analyzers.host_compatibility import ExpressionHostScore, ExpressionCompatibilityReport  # noqa: F401
except Exception:
    pass

# ── P6: Druggability Report ────────────────────────────────────────────────────
try:
    from analyzers.druggability_analyzer import BindingPocket, DruggabilityReport  # noqa: F401
except Exception:
    pass

# ── Ph2: Conformational Analysis ─────────────────────────────────────────────
try:
    from analyzers.conformational_analyzer import FlexibilityMode, ConformationalAnalysis  # noqa: F401
except Exception:
    pass

# ── P9: Inverse Folding ───────────────────────────────────────────────────────
try:
    from analyzers.inverse_folder import SynthesisFlag, InverseDesign  # noqa: F401
except Exception:
    pass

# ── P10: Therapeutic Decision Simulator ──────────────────────────────────────
try:
    from analyzers.therapeutic_simulator import ModalityScore, TherapeuticDecisionReport  # noqa: F401
except Exception:
    pass

# ── P11: Antigen Discovery Mode ───────────────────────────────────────────────
try:
    from analyzers.antigen_discovery_analyzer import AntigenCandidate, AntigenDiscoveryReport  # noqa: F401
except Exception:
    pass

# ── Layer 1: Allosteric Network ───────────────────────────────────────────────
try:
    from analyzers.allosteric_analyzer import AllostericHub, AllostericNetwork  # noqa: F401
except Exception:
    pass

# ── Layer 1: Epistasis Analysis ───────────────────────────────────────────────
try:
    from analyzers.epistasis_analyzer import EpistasisPair, EpistasisReport  # noqa: F401
except Exception:
    pass

# ── Layer 1: Binding Energy Estimation ───────────────────────────────────────
try:
    from analyzers.binding_energy_estimator import BindingEnergyEstimate  # noqa: F401
except Exception:
    pass

# ── Layer 1: Observable Predictor ────────────────────────────────────────────
try:
    from analyzers.observable_predictor import ObservablesPrediction  # noqa: F401
except Exception:
    pass

# ── Layer 2: Active Learning Advisor ─────────────────────────────────────────
try:
    from analyzers.active_learning_advisor import MutationRecommendation  # noqa: F401
except Exception:
    pass

# ── Layer 3: Glycan Analysis ─────────────────────────────────────────────────
try:
    from analyzers.glycan_analyzer import GlycanSite, GlycanAnalysis  # noqa: F401
except Exception:
    pass

# ── Layer 3: ADMET Prediction ─────────────────────────────────────────────────
try:
    from analyzers.admet_analyzer import ADMETReport  # noqa: F401
except Exception:
    pass

# ── Layer 3: Directed Evolution ───────────────────────────────────────────────
try:
    from analyzers.directed_evolution_analyzer import DirectedEvolutionPlan  # noqa: F401
except Exception:
    pass

# ── Layer 3: PROTAC Feasibility ───────────────────────────────────────────────
try:
    from analyzers.protac_analyzer import PRORTACFeasibilityReport  # noqa: F401
except Exception:
    pass

# ── Layer 3: Selectivity Landscape ───────────────────────────────────────────
try:
    from analyzers.selectivity_analyzer import SelectivityLandscape  # noqa: F401
except Exception:
    pass

# ── Layer 3: Fragment Hotspot Map ─────────────────────────────────────────────
try:
    from analyzers.fragment_hotspot_analyzer import FragmentHotspotMap  # noqa: F401
except Exception:
    pass

# ── Layer 3: AAV Gene Therapy Designer ───────────────────────────────────────
try:
    from analyzers.aav_designer import AAVDesignReport  # noqa: F401
except Exception:
    pass

# ── Layer 3: CRISPR Integration Planner ──────────────────────────────────────
try:
    from analyzers.crispr_designer import CRISPRIntegrationPlan  # noqa: F401
except Exception:
    pass

# ── Layer 3: Genetic Circuit Architect ───────────────────────────────────────
try:
    from analyzers.genetic_circuit_architect import GeneticCircuitDesign  # noqa: F401
except Exception:
    pass

# ── Layer 3+: Covalent Inhibitor Designer ────────────────────────────────────
try:
    from analyzers.covalent_inhibitor_designer import WarheadOption, CovalentSite, CovalentInhibitorDesign  # noqa: F401
except Exception:
    pass

# ── Layer 3+: Protein Engineering Advisor ────────────────────────────────────
try:
    from analyzers.engineering_strategy_recommender import StabilizationStrategy, EngineeringStrategyReport  # noqa: F401
except Exception:
    pass

# ── Layer 3+: H-Bond Network Analyzer ────────────────────────────────────────
try:
    from analyzers.hbond_network_analyzer import HBond, HBondCluster, HBondNetwork  # noqa: F401
except Exception:
    pass

# ── V3: ETP Mapper (Beratan-Onuchic) ─────────────────────────────────────────
try:
    from analyzers.etp_mapper import ETPResidue, ETPEdge, ETPAnalysis  # noqa: F401
except Exception:
    pass

# Resolve forward references
ProteinRecord.model_rebuild()
