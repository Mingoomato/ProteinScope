# ProteinScope — Biological Consequence Engine

> *"Don't just ask what a protein is. Ask what happens when it changes."*

**ProteinScope** is a full-stack AI-powered protein analysis platform that transforms a single gene symbol or UniProt accession into a comprehensive, multi-layer scientific report — integrating 20+ live biological databases, 46 specialized analysis modules, ESM-2 deep learning, AlphaFold structural analysis, Gemini 2.5 Pro AI synthesis, and an interactive 3D molecular viewer.

Built by a single developer as a solo engineering project, ProteinScope evolved from a simple protein lookup tool into what the developer calls a **Biological Consequence Engine**: a platform that answers not just *"what is this protein?"* but *"what will happen if I mutate it, drug it, edit it, or engineer it?"*

---

## Table of Contents

1. [Why This Exists](#why-this-exists)
2. [What ProteinScope Does](#what-proteinscope-does)
3. [The 4-Layer Architecture](#the-4-layer-architecture)
4. [Complete Feature Inventory](#complete-feature-inventory)
5. [The 3-AI Brainstorming Session](#the-3-ai-brainstorming-session)
6. [Technical Architecture](#technical-architecture)
7. [Database Integrations](#database-integrations)
8. [Development Journey & Challenges](#development-journey--challenges)
9. [Getting Started](#getting-started)
10. [Project Statistics](#project-statistics)
11. [Research Foundation](#research-foundation)
12. [Roadmap](#roadmap)

---

## Why This Exists

The modern life sciences researcher faces a paradox: there is more biological data available than at any point in history, yet transforming that data into actionable insight still requires jumping between dozens of tools, databases, and file formats — losing context at every step.

A typical protein investigation workflow looks like this:

```
UniProt → download FASTA
→ BLAST for homologs
→ AlphaFold for structure
→ PharmGKB for drug interactions
→ ClinVar for clinical variants
→ STRING for PPI network
→ PhosphoSitePlus for PTMs
→ ...manually integrate everything
→ write it up for a meeting
```

Each step involves a context switch. Each tool has a different interface, a different data format, a different update cadence. By the time you've assembled the picture, a colleague has moved on, the meeting has happened, or the drug target opportunity has been scooped.

**ProteinScope collapses this entire workflow into a single query.**

Type a gene name. Get a structured, cited, AI-synthesized report covering molecular function, structural analysis, evolutionary conservation, drug interactions, clinical variants, PTM logic, protein engineering options, CRISPR/AAV design parameters, covalent inhibitor sites, and a Gemini 2.5 Pro narrative — all in under 30 seconds.

---

## What ProteinScope Does

### Core Query Flow

```
User types: "EGFR" or "P00533"
                    │
                    ▼
         Intent Router (8 query types)
                    │
         ┌──────────┼──────────────┐
         ▼          ▼              ▼
    UniProt      AlphaFold     20+ DB fetchers
    (primary)    (structure)   (parallel async)
         │          │              │
         └──────────┴──────────────┘
                    │
                    ▼
           46 Analysis Modules
           (sequential + parallel)
                    │
                    ▼
         Gemini 2.5 Pro Synthesis
                    │
                    ▼
        Real-time SSE streaming to UI
        Mol* 3D viewer | Cytoscape PPI
        Downloadable PDF/Markdown report
```

### Query Modes

| Mode | Trigger | Output |
|------|---------|--------|
| **Standard** | Gene name / accession | Full protein report |
| **Concept Search** | "proteins involved in X" | Semantic FAISS search over indexed proteins |
| **RNAi Design** | "knockdown BRCA1" | siRNA/shRNA sequences, off-target analysis |
| **Reverse Genetics** | "what phenotype if I knock out X?" | Pathway impact, model organism context |
| **Compatibility** | "express KRAS in yeast?" | Host compatibility scoring, codon optimisation |
| **Antigen Discovery** | Tumor type + modality | CAR-T / BiTE / ADC antigen candidates from HPA |
| **Batch Mode** | Multiple genes | Sequential multi-protein reports |
| **Pathway Planning** | Disease/pathway context | Intervention hierarchy |

---

## The 4-Layer Architecture

The platform is organized around a **Biological Consequence Engine** — four layers that progressively transform molecular facts into therapeutic action.

```
┌─────────────────────────────────────────────────────────────────┐
│  LAYER 4 — Therapeutic Strategy                                 │
│  Modality scoring (8 drug classes) · Clinical trial context    │
│  Antigen discovery · CAR-T / BiTE / ADC / PROTAC feasibility   │
├─────────────────────────────────────────────────────────────────┤
│  LAYER 3 — Engineering Feasibility                              │
│  Glycan sites · ADMET · Directed evolution · CRISPR guides     │
│  AAV vector design · Genetic circuits · Fragment hotspots      │
│  Covalent inhibitor design · Protein thermostabilisation       │
├─────────────────────────────────────────────────────────────────┤
│  LAYER 2 — Predictive Mechanism                                 │
│  Epistasis (DMS interpretation) · Binding energy (MM/GBSA)    │
│  Observable prediction (NMR/FRET/HDX) · Active learning        │
├─────────────────────────────────────────────────────────────────┤
│  LAYER 1 — Molecular State Model                                │
│  Allosteric network (betweenness centrality) · PTM logic FSM   │
│  IDP/LLPS (FuzDrop/IUPred3) · Metalloenzyme coordination      │
│  H-bond network · Conformational ANM/GNM · Inverse folding     │
└─────────────────────────────────────────────────────────────────┘
```

This architecture was designed in a **3-way AI brainstorming session** (see [below](#the-3-ai-brainstorming-session)).

---

## Complete Feature Inventory

### Layer 1 — Molecular State Model

| Module | Description | Key Citations |
|--------|-------------|---------------|
| **ESM-2 Variant Scoring** | ΔLog-likelihood for every clinical variant; D3 heatmap visualization | Rives 2021 PNAS |
| **PTM Logic Engine** | State-transition graph of phosphorylation/ubiquitination dependencies | PhosphoSitePlus; UniProt |
| **IDP / LLPS Analysis** | Intrinsically disordered regions; liquid-liquid phase separation propensity | IUPred3; FuzDrop (Erdős 2021) |
| **Metalloenzyme Analysis** | Metal coordination geometry; pKa shifts (PROPKA); cofactor binding | Andreini 2008; Jumper 2021 |
| **Conformational Analysis** | ANM/GNM normal mode analysis (ProDy); collective motion modes | Yang 2008 |
| **Allosteric Network** | Cα contact graph → betweenness centrality; PhaSepDB LLPS overlay | Dokholyan 2002; You 2020 |
| **H-Bond Network** | AlphaFold structure → full H-bond network; buried cluster detection | McDonald 1994; Baker 1984 |
| **Inverse Folding** | ESM-IF1 sequence design from AlphaFold structure; 3-gate filtering | Hsu 2022; Kozak 1986 |

### Layer 2 — Predictive Mechanism

| Module | Description | Key Citations |
|--------|-------------|---------------|
| **Epistasis Analysis** | Pairwise epistasis index from DMS data; MaveDB integration | Domingo 2019 |
| **Binding Energy Estimation** | MM/GBSA-proxy ΔG for all drug interactions; pocket volume scaling | Gilson 2007; Halgren 2009 |
| **Observable Prediction** | NMR suitability (MW/TROSY), FRET Cys pair candidates, HDX fast-exchange | Kay 2011; Stryer 1978 |
| **Active Learning Advisor** | UCB acquisition for next mutation experiment; MaveDB overlap | Frazier 2018 |

### Layer 3 — Engineering Feasibility

| Module | Description | Key Citations |
|--------|-------------|---------------|
| **Glycan Analysis** | N/O-glycosylation sequon prediction; GlyConnect evidence | Apweiler 1999; Gupta 2004 |
| **ADMET Prediction** | Lipinski rule-of-5, BBB penetration, hepatotoxicity flags | Lipinski 2001; Egan 2000 |
| **Directed Evolution** | DMS-guided library design; error-prone PCR / site-saturation strategies | Romero 2009 |
| **PROTAC Feasibility** | E3 ligase compatibility; linker chemistry; degradation score | Sakamoto 2001; Bekes 2022 |
| **Selectivity Landscape** | Paralog pocket-identity scoring; cross-reactivity risk tiers | Hopkins 2002 |
| **Fragment Hotspot Map** | Fragment-based druggability; site-specific buriedness scoring | Erlanson 2016 |
| **AAV Gene Therapy Designer** | Serotype selection (AAV1–AAV9); payload capacity check; promoter | Grieger 2005 |
| **CRISPR Integration Planner** | Guide RNA design (SpCas9/Cas12); on-target scoring; base editing | Doench 2016; Anzalone 2019 |
| **Genetic Circuit Architect** | iGEM BioBrick topology selection; promoter/RBS/backbone pairing | Brophy 2014; Gardner 2000 |
| **Covalent Inhibitor Designer** | Nucleophilic hotspot scanning (Cys/Ser/Lys/Tyr/His); warhead chemistry | Bauer 2015; Baillie 2016 |
| **Protein Engineering Advisor** | ΔTm estimation; Pro/disulfide/buried-polar strategies; expression tips | Guerois 2002; Eijsink 2004 |

### Layer 4 — Therapeutic Strategy

| Module | Description | Key Citations |
|--------|-------------|---------------|
| **Therapeutic Simulator** | 8-modality scoring (small molecule, biologic, PROTAC, gene therapy, etc.) | Pammolli 2011; Hay 2014 |
| **Antigen Discovery** | Tumor vs. normal HPA expression; CAR-T/BiTE/ADC scoring | Carter 2001; Klebanoff 2016 |

### Infrastructure Modules

| Module | Description |
|--------|-------------|
| **Sequence Aligner** | BLOSUM62 + NUC44 pairwise; MSA; conservation scores |
| **Phylogenetic Analyzer** | Cross-species conservation; ortholog alignment |
| **RNAi Analyzer** | siRNA/shRNA design (Tuschl rules, Reynolds scoring, off-target) |
| **Compatibility Analyzer** | Cross-organism expression compatibility; TM-score comparison |
| **Reverse Genetics** | Phenotype prediction from knockout; model organism mapping |
| **Embedding Search** | FAISS-indexed protein concept search (sentence-transformers) |
| **Literature RAG** | PubMed retrieval-augmented generation; citation-grounded answers |
| **Pathway Integrator** | KEGG + Reactome + MetaCyc integration |
| **PGx Reporter** | PharmGKB evidence-graded pharmacogenomics variants |
| **Druggability Analyzer** | fpocket pocket detection; ChEMBL cross-reference |
| **Domain Scanner** | HMMER3 Pfam scan; domain architecture |
| **Structure Comparator** | TM-score; RMSD; structural homology |
| **Session Manager** | SQLite-backed session context; project board sidebar |

---

## The 3-AI Brainstorming Session

On **2026-03-19**, an unconventional design session took place: three AI systems were convened as domain specialists to architect the next phase of ProteinScope's development.

### Setup & Methodology

Rather than using AI as a code generator, the session used AI in **expert advisory roles** — each system was assigned two specialist personas and asked to submit a written opening position paper before the discussion began.

```
┌─────────────────────────────────────────────────────────────────┐
│           TRIPARTITE AI BRAINSTORMING SESSION                   │
│                    2026-03-19                                   │
│                                                                 │
│  Google Gemini 3        OpenAI GPT-5.4 (Codex)    Claude 4.6  │
│  ─────────────────       ──────────────────────    ──────────  │
│  Biochemist            Computational               Biotech     │
│  Quantum Biologist     Bioinformatician            Medicinal   │
│                        Biophysicist                Chemist     │
│                                                                 │
│  Facilitator & Scribe: Claude 4.6                              │
└─────────────────────────────────────────────────────────────────┘
```

### How the Session Ran

**Phase 1 — Opening Position Papers**

Each AI was asked: *"From your specialist perspective, what is ProteinScope's single most significant scientific blind spot?"*

- **Gemini** (Biochemist): *"ProteinScope treats every protein as a single static statue. Proteins are dynamic statistical ensembles. The platform is giving researchers a photo of a moving dancer, frozen mid-step."*

- **Codex** (Bioinformatician): *"The platform has no mechanism for representing epistasis — the non-additive interactions between mutations. Every variant is scored in isolation, which is scientifically incorrect for proteins with cooperative binding or allosteric regulation."*

- **Claude** (Medicinal Chemist): *"The platform can tell you everything about a protein, but it cannot tell a drug discovery team which of their next ten experiments will give the most information. There is no experimental design layer."*

**Phase 2 — Cross-Examination**

Each AI challenged the others' proposals with scientific rigor:

- Codex challenged Gemini's quantum tunneling proposal: *"Proton-coupled electron transfer probability maps require quantum chemistry calculations. What is the fallback when no crystal structure exists? A web service cannot run DFT."* → Resolution: implement as a gene-symbol lookup for known PCET enzymes with inline citations; flag rather than calculate.

- Gemini challenged Claude's CRISPR design feature: *"Guide RNA design without structural context of the target site ignores chromatin accessibility. You are designing in sequence space while ignoring 3D context."* → Resolution: combine CRISPR guide scoring with AlphaFold domain boundaries; note chromatin data as out of scope.

- Claude challenged Gemini's H-bond network proposal: *"Which H-bonds matter for drug binding vs. fold stability? Without that distinction, you're showing researchers noise."* → Resolution: cluster by proximity and classify as backbone-backbone vs. sidechain-sidechain; flag buried clusters specifically.

**Phase 3 — Convergence on Principles**

Three principles emerged from the disagreements:

1. **Every number must show its evidence grade and source.** The Feature Confidence Framework (`EvidenceGrade` + `DataProvenance`) became the architectural prerequisite for all new features.

2. **Actionable output over completeness.** Every layer must end with something a researcher can *do* — a sequence to synthesize, a guide RNA to order, a warhead to try — not just a number to look at.

3. **Communicate uncertainty as clearly as confidence.** The quantum tunneling flag is not a failure of scope; it is the platform being honest about what physics requires. Showing a "flag" rather than a fake probability is the scientifically correct choice.

**Phase 4 — Architecture Proposal**

The 4-layer Biological Consequence Engine architecture emerged from this negotiation — not designed top-down, but negotiated bottom-up from three different scientific worldviews:

| Layer | Origin |
|-------|--------|
| Layer 1: Molecular State Model | Gemini's dynamic protein ensemble argument |
| Layer 2: Predictive Mechanism | Codex's epistasis + biophysics argument |
| Layer 3: Engineering Feasibility | Claude's experimental design argument |
| Layer 4: Therapeutic Strategy | Consensus across all three |

The full meeting transcript is preserved in `TRIPARTITE_BRAINSTORM_MINUTES.md`.

### Why This Methodology Matters

This session demonstrated that AI systems can be orchestrated as **domain-specific advisory panels** — a methodology with implications beyond software development.

The key insight: **asking "what is wrong with our approach?" produces more rigorous design than asking "what should we build?"** Disagreements between systems force justification of every design decision. Features that survive cross-examination are features that will survive peer review.

This is now the standard design process for ProteinScope: major architectural decisions require a multi-AI review session before implementation begins.

---

## Technical Architecture

### Backend Stack

```
┌──────────────────────────────────────────────────────────────┐
│  FastAPI (async)  ·  Uvicorn ASGI  ·  Server-Sent Events    │
├──────────────────────────────────────────────────────────────┤
│  Python 3.11+                                                │
│  Pydantic v2 (64-field ProteinRecord + 100+ sub-models)     │
│  asyncio + httpx (parallel DB fetching, graceful degrade)   │
├──────────────────────────────────────────────────────────────┤
│  Google Gemini 2.5 Pro (AI synthesis & narrative)           │
│  ESM-2 650M (Meta AI, mutation fitness scoring via ΔLL)     │
│  ESM-IF1 (Meta AI, inverse folding / sequence design)       │
│  sentence-transformers (semantic protein concept search)    │
│  FAISS (vector similarity search, protein embedding index)  │
├──────────────────────────────────────────────────────────────┤
│  ProDy (ANM/GNM normal mode analysis)                       │
│  NetworkX (allosteric betweenness centrality graphs)        │
│  Biopython · PyHMMER · biotite (bioinformatics core)        │
│  ReportLab (PDF report generation)                          │
├──────────────────────────────────────────────────────────────┤
│  SQLite (session manager, 7-day analysis cache)             │
│  Docker (CPU image + CUDA GPU image)                        │
└──────────────────────────────────────────────────────────────┘
```

### Frontend Stack

```
┌──────────────────────────────────────────────────────────────┐
│  Vanilla JS (ES2022)  ·  Tailwind CSS (JIT via CDN)         │
│  Server-Sent Events (real-time streaming progress)          │
├──────────────────────────────────────────────────────────────┤
│  Mol* Viewer (3D molecular visualization, RCSB-standard)    │
│  Cytoscape.js + CoSE-Bilkent (PPI network graph)            │
│  D3.js (ESM-2 fitness landscape heatmaps)                   │
├──────────────────────────────────────────────────────────────┤
│  Single-page, zero-framework  ·  ~5,000 line JS/HTML        │
│  17 specialized report display functions                    │
│  Evidence badge system (experimental/computational/AI)      │
│  PDF + Markdown export                                      │
└──────────────────────────────────────────────────────────────┘
```

### The Evidence Framework

Every computed value in ProteinScope carries a **DataProvenance** object. This is the platform's core scientific commitment: *uncertainty is first-class data*.

```python
class EvidenceGrade(str, Enum):
    EXPERIMENTAL = "experimental"    # Crystal structure, wet-lab assay, ClinVar
    COMPUTATIONAL = "computational"  # AlphaFold, ESM-2, PROPKA
    LITERATURE    = "literature"     # Extracted from PubMed
    AI_GENERATED  = "ai_generated"   # Gemini synthesis

class DataProvenance(BaseModel):
    source: str                           # "UniProt REST v2", "ESM-2 650M"
    evidence_grade: EvidenceGrade
    confidence_interval: Optional[str]    # "±2 kcal/mol (MM/GBSA proxy)"
    scientific_caveat: Optional[str]      # "AlphaFold; no experimental structure"
    method: Optional[str]
```

In the UI, every data point renders a colored evidence badge:
- 🟢 **Experimental** — from crystal structures, wet-lab assays, clinical databases
- 🔵 **Computational** — from AlphaFold, ESM-2, or structure prediction
- 🟡 **Literature** — extracted from PubMed papers
- 🟣 **AI-generated** — Gemini 2.5 Pro narrative synthesis

### The Analyzer Pattern

All 46 analysis modules follow a strict interface contract:

```python
async def run_X_analysis(
    gene: str,
    # ...domain-specific inputs...
    step_cb=None,            # Optional SSE progress callback
) -> Optional[XReport]:     # Never raises; returns None on failure
    try:
        # Step 1: fetch/compute
        # Step 2: ...
        # Step N: Gemini synthesis
        return XReport(...)
    except Exception:
        return None           # Graceful degradation — platform never crashes
```

This pattern ensures that adding a new analysis module never risks breaking existing functionality. A broken analyzer returns `None`; the ProteinRecord field is simply empty. The report renders without that section.

---

## Database Integrations

ProteinScope queries **20+ biological databases** in parallel on every protein lookup:

| Database | Data Retrieved | Notes |
|----------|---------------|-------|
| **UniProt REST v2** | Sequence, function, domains, features, PTMs, subcellular location | Required |
| **AlphaFold EBI v4** | PDB structure, per-residue pLDDT confidence | Graceful skip |
| **NCBI ClinVar** | Clinical variants, pathogenicity classifications | Graceful skip |
| **NCBI Gene** | Gene ID, gene summary, aliases | Graceful skip |
| **PDB (RCSB)** | Experimental structure IDs | Graceful skip |
| **PharmGKB** | Pharmacogenomics variants, drug-gene evidence levels | Graceful skip |
| **DGIdb v5** | Drug-gene interactions | Graceful skip |
| **ChEMBL** | Drug molecules, SMILES, mechanism of action | Graceful skip |
| **STRING v11** | Protein-protein interaction network | Graceful skip |
| **KEGG** | Metabolic/signaling pathways | Graceful skip |
| **Reactome** | Signaling pathways, disease annotations | Graceful skip |
| **MetaCyc** | Enzymatic reactions, metabolic context | Graceful skip |
| **PhosphoSitePlus** | Experimentally validated PTM sites | Optional API key |
| **IUPred3** | Per-residue disorder prediction | Graceful skip |
| **FuzDrop** | LLPS propensity, droplet-forming regions | Graceful skip |
| **Human Protein Atlas** | Tissue/tumor/subcellular expression | Graceful skip |
| **ClinicalTrials.gov** | Active clinical trials for the target | Graceful skip |
| **PubMed** | Citation metadata | Graceful skip |
| **MaveDB** | Deep mutational scanning datasets | Graceful skip |
| **PhaSepDB** | Phase separation; condensate partners | Graceful skip |
| **GlyConnect** | Glycoprotein database; glycan evidence | Graceful skip |

All fetchers follow the same contract: **return empty on failure, never propagate exceptions**.

---

## Development Journey & Challenges

### The Origin

ProteinScope began as a tool to answer a specific personal question: *"If I have a candidate drug target, what do I actually need to know about it before committing three months to a project?"*

The first version was a Python script that called UniProt and returned a text summary. It grew.

### Challenge 1: Scope Creep as a Structural Problem

Every feature added raised three new scientific questions. Adding drug interactions surfaced the need for pharmacogenomics context. Adding pharmacogenomics surfaced the need for variant scoring. Adding variant scoring surfaced the need for epistasis analysis. The architecture had to be designed to accommodate features whose names didn't exist when the project started.

**Solution**: The **Analyzer Pattern** — a strict interface (`async def run_X(..., step_cb=None) -> XReport | None`, never raises, always returns, Gemini synthesis as the final step). This meant each new capability slotted in without touching existing code. Adding a new analyzer requires touching exactly four files: the analyzer file itself, `models.py`, `query_engine.py`, and `index.html`. Nothing else needs to change.

### Challenge 2: 20 APIs, 20 Failure Modes

Running 20+ async database queries simultaneously against public research APIs is operationally complex. Each API has:
- Different rate limits (some per-IP, some per-key, some undocumented)
- Different response schemas (JSON, XML, TSV, binary formats)
- Different error conventions (many APIs return HTTP 200 with an error body)
- Different uptime characteristics (planned maintenance windows, aggressive crawl limits)

During development, queries that worked consistently for weeks would suddenly start returning 503s from one database while all others continued normally.

**Solution**: A **7-day SQLite cache** and strict **graceful degradation**. Every fetcher catches all exceptions and returns an empty structure. The UI never shows a 500 error; it shows a partial report with exactly the data that was available. A researcher gets 90% of the information even when three databases are simultaneously unavailable.

### Challenge 3: The Citation Discipline

Every numeric threshold in the codebase must be accompanied by an inline citation to a peer-reviewed paper. This was non-negotiable from early in the project.

```python
# 8 Å contact threshold for allosteric network construction
# Citation: Dokholyan NV et al. 2002 J Mol Biol doi:10.1006/jmbi.2001.5332
_CONTACT_CUTOFF_ANGSTROM = 8.0

# pLDDT < 70 → locally disordered region
# Citation: Jumper J et al. 2021 Nature doi:10.1038/s41586-021-03819-2
_LOW_PLDDT_THRESHOLD = 70.0

# Epistasis index threshold for "high-risk" combination
# Citation: Domingo J et al. 2019 Nat Commun doi:10.1038/s41467-019-09178-5
_HIGH_RISK_INDEX = 3.0
```

This discipline forced a rigorous examination of every number in the codebase. Several initially "obvious" thresholds were changed after the literature review revealed the actual experimentally-validated values. It also meant that every module is, in effect, a literature review as much as a code module.

### Challenge 4: Pydantic v2 Forward References at Scale

The ProteinRecord model grew from a handful of fields to **64 fields** spanning seven nested model types across multiple files. Pydantic v2's forward-reference resolution requires all referenced types to be registered before `model_rebuild()` is called.

**Solution**: A specific pattern for adding new models — try/except import blocks that allow graceful degradation when analyzer modules are incomplete:

```python
# Pattern established after several production failures:
# try/except import before model_rebuild() allows graceful degradation
# when analyzer modules are incomplete or import fails
try:
    from analyzers.new_module import NewModel
except Exception:
    pass  # model_rebuild() uses string annotations as fallback

ProteinRecord.model_rebuild()  # MUST remain the last line of models.py
```

Breaking this pattern (e.g., placing imports after `model_rebuild()`) caused silent failures where the ProteinRecord silently lost fields. This rule is now the most-documented constraint in the entire codebase.

### Challenge 5: Real-Time Streaming Architecture

Researchers expect to see progress. A 30-second analysis with no feedback feels like a crash. Designing a streaming progress architecture that was informative but not overwhelming required significant iteration.

**Solution**: A `step_cb` callback pattern. Every analysis function accepts an optional async callback. The analysis modules report meaningful steps ("downloading structure", "detecting h-bonds", "gemini synthesis") and the infrastructure above decides whether to actually stream them. This kept the analysis modules completely independent of the transport layer — the same module can run in a web server (SSE), a CLI (print), or a batch job (logging) with no code changes.

### Challenge 6: AI-Assisted Development at Scale

The development process itself used AI assistance (Claude) to accelerate implementation of the 46 analysis modules. During the session that implemented Layers 1–3, agent sub-processes for groups A3, A4, and A5 hit concurrent rate limits and returned failure messages instead of code.

Rather than treating this as a blocker, the session documented which agents completed and which failed, preserved the execution state to a persistent memory system, and the primary session implemented the remaining modules directly. **The development toolchain was itself subject to the same reliability engineering principles as the platform being built.**

This experience reinforced the value of the graceful degradation philosophy: a development process that crashes on partial failure is as problematic as a production system that crashes on partial data.

### Challenge 7: Context Management Across Long Sessions

As the codebase grew beyond 50,000 lines, maintaining development context across sessions became a fundamental challenge. AI-assisted development sessions have context windows — when a session ran out of context mid-execution, the work could not continue without reconstructing the state.

**Solution**: A **persistent memory system** — structured markdown files categorizing memories by type (user context, feedback, project state, references), with a MEMORY.md index that loads automatically into each new session. When a session was compressed mid-execution, the next session read the memory files, verified which files existed on disk, and continued from exactly where the previous session left off.

This was not a solution designed upfront. It emerged from necessity and was iteratively refined. The memory system now documents over 15 distinct project facts, including architectural invariants, citation requirements, and specific failure modes to avoid.

### The Philosophical Challenge

The deeper challenge was a scientific one: **where does a tool end and a research platform begin?**

Tools answer specific questions. Research platforms help researchers ask better questions. ProteinScope started as a tool. The 3-AI brainstorming session was the moment it became a research platform — when the design process itself required the kind of multi-disciplinary synthesis that the platform aims to enable.

The quantum tunneling flag is a small example of this distinction. It doesn't calculate tunneling probabilities (that requires quantum chemistry beyond what's feasible in a web service). Instead, it flags proteins where proton-coupled electron transfer has been experimentally studied, with a citation to Hammes-Schiffer 2006. **The platform communicates what it cannot do as clearly as what it can do.** That design decision was itself a research outcome of the 3-AI session — not a compromise, but a principled position about the difference between scientific honesty and feature completeness.

---

## Getting Started

### Prerequisites

- Python 3.11+
- Docker (recommended) or virtual environment
- Google Gemini API key

### Quick Start (Docker — recommended)

```bash
git clone https://github.com/your-username/protein-finder-project
cd protein-finder-project

# Copy environment template
cp proteinscope/.env.example proteinscope/.env
# Edit .env and add your GEMINI_API_KEY

# CPU build
docker-compose up

# GPU build (requires NVIDIA container toolkit)
docker-compose -f docker-compose.gpu.yml up

# Open http://localhost:8000
```

### Local Development

```bash
cd proteinscope
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt

echo "GEMINI_API_KEY=your_key_here" > .env
python main.py
```

### Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `GEMINI_API_KEY` | **Yes** | Google Gemini 2.5 Pro API key |
| `PHOSPHOSITE_API_KEY` | No | PhosphoSitePlus (enhances PTM data) |
| `BRENDA_EMAIL` | No | BRENDA enzyme database |
| `OMIM_API_KEY` | No | OMIM disease database |
| `PHARMGKB_API_KEY` | No | PharmGKB (falls back to public tier) |

### Example Queries

```
EGFR          → Lung cancer kinase: full oncology context, TKI interactions, CRISPR guides
TP53          → Tumour suppressor: variant landscape, PROTAC feasibility, therapeutic simulation
PCSK9         → Cholesterol target: druggability, covalent inhibitor sites, clinical trials
APP           → Alzheimer's amyloid: LLPS analysis, condensate partners, AAV design
BRCA1         → DNA repair: epistasis landscape, directed evolution, active learning advice
BTK           → B-cell kinase: ibrutinib context, warhead chemistry, selectivity risks
```

---

## Project Statistics

| Metric | Value |
|--------|-------|
| **Lines of Python (analyzers + core)** | ~21,000 |
| **Lines of JavaScript / HTML** | ~5,000 |
| **Analysis Modules** | 46 |
| **Async Database Fetchers** | 28 |
| **ProteinRecord Fields** | 64 |
| **External Databases Integrated** | 20+ |
| **Peer-reviewed citations in code** | 80+ |
| **Supported query modes** | 8 |
| **Report sections per query** | 35+ |
| **AI models integrated** | 3 (Gemini 2.5 Pro, ESM-2 650M, ESM-IF1) |
| **AI systems in architecture session** | 3 (Gemini 3, Codex GPT-5.4, Claude 4.6) |
| **Docker images** | 2 (CPU + CUDA GPU) |
| **Development model** | Solo, continuous iteration |

---

## Research Foundation

Every threshold, formula, and interpretation in the codebase carries an inline citation. Selected key papers:

**Structural Biology**
- Jumper J et al. 2021 *AlphaFold2*. Nature 596:583. doi:10.1038/s41586-021-03819-2
- Hsu C et al. 2022 *ESM-IF1 inverse folding*. ICML. doi:10.1101/2022.04.10.487779
- McDonald IK & Thornton JM 1994 *HBPLUS*. J Mol Biol 238:777. doi:10.1006/jmbi.1994.1334
- Dokholyan NV et al. 2002 *Protein folding networks*. J Mol Biol 317:855. doi:10.1006/jmbi.2001.5332

**Protein Language Models**
- Rives A et al. 2021 *ESM-1*. PNAS 118:e2016239118. doi:10.1073/pnas.2016239118
- Lin Z et al. 2023 *ESM-2 scaling*. Science 379:1123. doi:10.1126/science.ade2574

**Drug Discovery**
- Bauer RA 2015 *Targeted covalent inhibitors*. Drug Discov Today 20:1061. doi:10.1016/j.drudis.2015.05.005
- Baillie TA 2016 *Covalent drug mechanisms*. Angew Chem 55:13408. doi:10.1002/anie.201601091
- Bekes M et al. 2022 *PROTACs in clinical trials*. Nat Rev Drug Discov 21:181. doi:10.1038/s41573-021-00371-6
- Pammolli F et al. 2011 *Drug R&D productivity*. Nat Rev Drug Discov 10:428. doi:10.1038/nrd3405

**Computational Biology**
- Domingo J et al. 2019 *Epistasis in protein evolution*. Nat Commun 10:1. doi:10.1038/s41467-019-09178-5
- Guerois R et al. 2002 *FoldX stability*. J Mol Biol 320:369. doi:10.1016/S0022-2836(02)00442-4
- Yang LW et al. 2008 *ProDy ANM*. Structure 16:449. doi:10.1016/j.str.2007.12.014

**Synthetic Biology**
- Gardner TS et al. 2000 *Genetic toggle switch*. Nature 403:339. doi:10.1038/35002131
- Brophy JAN & Voigt CA 2014 *Genetic circuit design*. Nat Methods 11:508. doi:10.1038/nmeth.2926
- Doench JG et al. 2016 *CRISPR guide optimisation*. Nat Biotechnol 34:184. doi:10.1038/nbt.3437

---

## Roadmap

### Next Session Agenda (deferred to next 3-AI brainstorming)

These features require architectural decisions about computational infrastructure (external binaries, ML training pipelines, QM calculations) and are deferred pending a dedicated design session:

1. **Full Quantum Tunneling Probability Map** — Bell/Wigner KIE predictions; decision on QM subprocess
2. **Radical Pair Spin-Dynamics Analyzer** — Cryptochrome family detection; scope decision
3. **Vibrational Mode Profiler** — MD trajectory requirements vs. static NMA
4. **Beratan-Onuchic ETP Mapper** — Electron tunneling pathways in redox enzymes
5. **Contrastive PPI Specificity Predictor** — Siamese encoder for cross-reactivity prediction
6. **Family-Specialized LoRA Models** — Fine-tuned ESM-2 for kinases/GPCRs/ion channels
7. **Hidden-State Conformational Classifier** — HMM/MSM state population analysis
8. **Full APBS/DelPhi Electrostatics Engine** — Surface potential maps; Docker APBS integration

### Infrastructure

- REST API for programmatic access (OpenAPI-documented)
- Batch pipeline mode (CSV input → bulk PDF reports)
- Multi-user collaborative project boards
- Custom embedding fine-tuning on user protein sets
- Integration with electronic lab notebooks (Benchling, LabArchives)
- WebSocket upgrade for lower-latency streaming

---

## Architecture Overview

```
                    ┌─────────────────────────────────┐
                    │         User Interface            │
                    │   Browser · SSE · Mol* · D3     │
                    └──────────────┬──────────────────┘
                                   │ HTTP / SSE
                    ┌──────────────▼──────────────────┐
                    │         FastAPI Backend           │
                    │   web/app.py · 15 endpoints      │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │       Query Engine               │
                    │   core/query_engine.py           │
                    │   Intent routing → 8 modes      │
                    └──────────────┬──────────────────┘
                     ┌─────────────┼─────────────┐
       ┌─────────────▼──┐  ┌───────▼──────┐  ┌──▼──────────────┐
       │  28 Fetchers   │  │ 46 Analyzers │  │  Gemini 2.5 Pro │
       │  (async, 7d    │  │  (4 layers,  │  │  (synthesis +   │
       │   SQLite cache)│  │   step_cb)   │  │   narrative)    │
       └────────────────┘  └──────────────┘  └─────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │       ProteinRecord              │
                    │   core/models.py · 64 fields    │
                    │   Pydantic v2 · EvidenceGrade   │
                    └─────────────────────────────────┘
```

---

## License

This project is a solo research and engineering portfolio work. Academic and research use is welcomed with attribution.

---

## Acknowledgements

**Biological Databases**: UniProt, AlphaFold EBI, NCBI (ClinVar, Gene), PharmGKB, DGIdb, ChEMBL, STRING, KEGG, Reactome, MetaCyc, PhosphoSitePlus, IUPred3, FuzDrop, Human Protein Atlas, ClinicalTrials.gov, MaveDB, PhaSepDB, GlyConnect.

**AI & ML Models**: ESM-2 and ESM-IF1 by Meta AI Research. Gemini 2.5 Pro by Google DeepMind. sentence-transformers by UKP Lab (Reimers & Gurevych 2019).

**Bioinformatics Libraries**: FastAPI, Pydantic, ProDy, Biopython, PyHMMER, biotite, NetworkX, FAISS, ReportLab.

**Visualization**: Mol* (RCSB), Cytoscape.js, D3.js.

**Architecture Session**: Google Gemini 3, OpenAI GPT-5.4 (Codex), Anthropic Claude 4.6.

---

*Built with obsessive attention to scientific rigor, one citation at a time.*
