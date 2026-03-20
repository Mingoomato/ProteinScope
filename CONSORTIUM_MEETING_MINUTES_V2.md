# ProteinScope Grand Consortium Meeting — Session II
## Multi-Departmental AI Brainstorming Conference
## Date: 2026-03-20 | Location: Virtual — ProteinScope Development Summit

---

## Participants

### Department A — Computational Bioinformatics & Biophysics (Codex Division)
| Name | Title | Background |
|---|---|---|
| **김박사 (Dr. Kim)** | Computational Bioinformatician | Stanford PhD Biology · MIT PhD Bioinformatics · Harvard PhD Information Engineering |
| **차박사 (Dr. Cha)** | Biophysicist | Oxford PhD Biophysics |
| **용박사 (Dr. Yong)** | Computational Biophysical Chemist | SNU PhD Bioinformatics · KAIST PhD Biophysical Informatics · KAIST PhD Computational Chemistry |

### Department B — Biochemistry & Quantum Biology (Gemini Division)
| Name | Title | Background |
|---|---|---|
| **노박사 (Dr. Noh)** | Biochemist | POSTECH PhD Biochemistry |
| **윤박사 (Dr. Yoon)** | Quantum Biologist | CALTECH PhD Quantum Mechanics · CALTECH PhD Quantum Biology |

### Department C — AI Biotechnology & Medicinal Chemistry (Claude Division)
| Name | Title | Background |
|---|---|---|
| **Claudin, PhD** | AI-Associated Biotech Developer | Mayo Clinic PhD — AI-Associated Biotechnology Development |
| **Dr. Konberg** | Medicinal Chemist | MIT Professor — Medicinal Chemistry |

---

## Chair's Opening Remarks

**Claudin:** Welcome to the second grand consortium meeting. Today's agenda is substantially larger than our first session. We now have a working platform — 46 analysis modules, 64-field ProteinRecord, live SSE streaming, Mol* 3D viewer — and we need to decide the next architectural leap. I'm going to open the floor immediately. The agenda has 12 items. We'll move fast.

The items fall into four clusters:
1. **Disease Association Layer** (Agenda #12 — new DB integrations, model design)
2. **Quantum Biology Track** (Agenda #1–4 — tunneling, radical pair, ETP, vibrational modes)
3. **ML Architecture Track** (Agenda #5–6 — contrastive PPI, LoRA family models)
4. **Engineering & UX** (Agenda #7–11 — conformational MSM, APBS, step skipping, UX bugs)

Let's begin.

---

## Session 1 — Disease Association Layer

**Chair: Claudin**

**Claudin:** Currently the platform shows KEGG H-numbers and OMIM stubs. For a protein like TLR2 we should be returning "associated with lepromatous leprosy (OMIM:246300), septic shock susceptibility, COVID-19 severity" — with evidence scores. We have six candidate databases. I want each of you to stake a position on which ones we use and why.

**김박사:** From a data integration standpoint, the answer is OpenTargets first, always. Their platform already aggregates genetic association (GWAS), somatic mutations (COSMIC), drug target evidence (ChEMBL), and literature co-mention — all pre-scored. Their GraphQL API is excellent. The `associatedDiseases` endpoint gives us everything in one call: disease name, EFO ontology ID, genetic association score, overall score. For a bioinformatics pipeline, you don't want to federate five databases when one federation layer already exists.

**용박사:** I agree on OpenTargets as the anchor, but I'd argue DisGeNET runs in parallel, not as a fallback. DisGeNET has 1.2 million gene-disease associations versus OpenTargets' more curated but narrower set. Their GDA (Gene-Disease Association) score is computed differently — it integrates variant-level, gene-level, and text-mining evidence with a Fuente score that's been validated against known Mendelian genes. For rare disease coverage specifically, DisGeNET outperforms OpenTargets. The Ensembl ID gap that Claudin flagged is real — we need to resolve Ensembl IDs. My recommendation: fetch from UniProt cross-references field (`xrefs → Ensembl`). It's already in the UniProt JSON we pull.

**차박사:** From a biophysical perspective, disease associations are most actionable when we can connect them to structural features. What I'd add to this layer is HPO (Human Phenotype Ontology) integration — not just as labels but as mechanistic connectors. HPO term HP:0001250 (Seizures) linked to an ion channel protein tells the platform something about the biophysical disease mechanism. The NCBI HPO API is free and returns HPO terms per gene with OMIM cross-references. This adds a phenotypic mechanistic layer that pure disease-name lists don't capture.

**노박사:** From biochemistry: what I'm missing from all these databases is the **molecular mechanism of pathogenicity**. OpenTargets gives you a score; it doesn't tell you *why* the association exists biochemically. For enzymatic targets, the association is usually loss-of-function (enzyme insufficiency) or gain-of-toxic-function. For receptor targets, it's signaling dysregulation. My recommendation is that after we pull the disease list, Gemini synthesizes a one-sentence mechanistic explanation per disease: "TLR2 gain-of-function → constitutive NF-κB → chronic inflammation in [disease]." That transforms a database lookup into actionable biochemical context.

**윤박사:** I want to flag something none of you have mentioned: **somatic vs germline distinction matters enormously** for quantum biology relevance. A germline variant in TLR2 is an innate immunity story. A somatic mutation in EGFR in NSCLC is a different story — it's about clonal evolution, selection pressure, and potentially altered electron transfer in the kinase domain (EGFR catalysis involves proton transfer steps). OpenTargets does distinguish somatic vs germline evidence scores. We should surface this distinction in the UI with different color coding — germline = amber, somatic = red, drug target = blue.

**Dr. Konberg:** As a medicinal chemist: the most clinically relevant database for my purposes is actually **ClinGen** (Clinical Genome Resource), which curates gene-disease validity classifications (Definitive, Strong, Moderate, Limited) based on experimental evidence. OpenTargets scores are good for discovery but ClinGen gives you the regulatory-grade confidence level that matters for therapeutic targeting decisions. If a gene has "Definitive" ClinGen classification for a disease, that's a much stronger therapeutic rationale than an OpenTargets score of 0.6. I'd add ClinGen as a fifth source, even though their API is less developer-friendly.

---

### Session 1 — Decisions

**Claudin (summarizing):**

| Decision | Resolution |
|---|---|
| Primary DB | OpenTargets (GraphQL) — anchor source, most comprehensive aggregated scoring |
| Secondary DB | DisGeNET (REST free tier) — rare disease coverage, independent GDA score |
| Tertiary | HPO API (NCBI) — phenotypic mechanistic layer |
| Fourth | ClinGen (Dr. Konberg's addition) — validity classification for therapeutic confidence |
| Canonical Ontology | **MONDO** — cross-ontology bridge; map EFO (OpenTargets) + DO (DisGeNET) + OMIM → MONDO IDs |
| Ensembl ID sourcing | UniProt cross-reference field (`xrefs → Ensembl`) — already fetched, zero new API calls |
| Evidence type coloring | germline = amber, somatic = red, drug target = blue, literature = gray |
| Gemini synthesis | 1-sentence mechanistic explanation per disease — no박사's recommendation |
| ClinGen validity | Add as badge: "Definitive / Strong / Moderate / Limited" on high-scoring associations |

**`DiseaseAssociation` model (finalized):**
```python
class DiseaseAssociation(BaseModel):
    disease_name: str
    mondo_id: Optional[str] = None          # MONDO cross-ontology ID
    efo_id: Optional[str] = None            # from OpenTargets
    omim_id: Optional[str] = None
    source_db: str                           # "OpenTargets" | "DisGeNET" | "ClinGen" | "HPO"
    overall_score: float                     # 0.0–1.0
    genetic_score: Optional[float] = None
    somatic_score: Optional[float] = None
    drug_score: Optional[float] = None
    evidence_type: str                       # "germline" | "somatic" | "drug_target" | "literature"
    clingen_validity: Optional[str] = None  # "Definitive" | "Strong" | "Moderate" | "Limited"
    hpo_terms: List[str] = Field(default_factory=list)
    mechanistic_summary: Optional[str] = None   # Gemini synthesis
    n_publications: Optional[int] = None
    provenance: Optional[DataProvenance] = None
```

**New files required:**
- `fetchers/opentargets.py` — GraphQL client
- `fetchers/disgenet.py` — REST client (free tier, rate limit: 1 req/sec)
- `fetchers/clingen.py` — ClinGen gene validity REST
- `analyzers/disease_association_analyzer.py` — deduplication + MONDO mapping + Gemini synthesis

---

## Session 2 — Quantum Biology Track

**Chair: 윤박사**

**윤박사:** I've been waiting for this session since the first meeting. The quantum biology features we deferred are not peripheral curiosities — they are increasingly recognized as central to enzyme catalysis, magnetoreception, and potentially neurotransmitter binding. Let me structure this: we have four deferred items. I'll frame each one with a concrete scientific question, then we debate feasibility.

---

### 2A — Quantum Tunneling Probability Map

**윤박사:** The central question: for a given enzyme with a known cofactor, what is the probability that the rate-limiting step involves quantum tunneling rather than classical over-barrier transfer? This is answerable computationally using Bell tunneling corrections (Bell 1935) and the Wigner correction factor (Wigner 1932). The KIE (kinetic isotope effect) — kH/kD — is the experimental observable; KIE > 7 at room temperature is strong evidence for tunneling. Can we predict this from sequence + structure alone?

**용박사:** Yes, with caveats. Here's the computational path: (1) identify the hydrogen donor-acceptor pair from the active site annotation. (2) Extract the donor-acceptor distance from AlphaFold PDB — must be < 3.5 Å for relevant tunneling probability. (3) Apply Bell tunneling correction: `Γ_Bell = (u/2) / sin(u/2)` where `u = hν‡/kT`, ν‡ is the imaginary frequency of the transition state. The problem: ν‡ requires QM/MM at minimum. However, for well-characterized enzyme families (alcohol dehydrogenase, DHFR, aromatic amine dehydrogenase), ν‡ values are tabulated in the literature. We can use a lookup table approach for the top 20 tunneling-documented enzymes and flag "tunneling likely" for structural homologs.

**차박사:** I want to be precise about what we can and cannot compute without MD. The Wigner correction gives order-of-magnitude estimates; it's not quantitative without the PES (potential energy surface). What we *can* do without MD: use the Swain-Schaad exponent as a tunneling indicator. If the enzyme belongs to a family where tunneling is documented (DHFR → tunneling at C4H of NADPH; Aromatic amine dehydrogenase → tryptamine oxidation), we flag it with the literature-supported KIE range, not a computed value. This is honest: "Tunneling documented in this enzyme class (KIE_exp = 8.2, Hay & Scrutton 2012)."

**김박사:** For the computational implementation, I'd use a two-tier architecture: Tier 1 — homology-based lookup (HMMER against tunneling-documented enzyme families, E-value < 1e-20 → flag with literature KIE); Tier 2 — structural geometry check (donor-acceptor distance from AlphaFold PDB using BioPython PDB parser, check for tunneling-relevant configuration). This requires: HMMER (already in requirements? No — add `hmmer` via subprocess), BioPython (installed), and a curated `data/tunneling_enzymes.json` lookup table of ~50 enzyme families.

**윤박사:** Accepted. Decision: **implement Tier 1 (HMMER lookup) + Tier 2 (geometry check) — flag with literature KIE, never fabricate a computed ν‡ value without QM/MM.**

---

### 2B — Radical Pair Spin-Dynamics Analyzer

**윤박사:** Cryptochrome-based magnetoreception (Ritz 2000, PNAS) involves radical pair formation in flavin cofactors. The spin-coherence lifetime determines the magnetic field sensitivity. For proteins that contain FAD/FMN cofactors and belong to the photolyase/cryptochrome family: can we flag them and estimate spin-coherence lifetime?

**용박사:** The spin Hamiltonian for a radical pair is: `H = Σ_i Σ_k a_ik S_i · I_k + Σ_i g_i μ_B B_0 · S_i` where a_ik is the hyperfine coupling constant and g_i is the g-factor. For flavin radicals (FADH•), the hyperfine coupling constants are well-characterized (Maeda 2012, Science). We can compute the expected coherence lifetime using the Schulten model (Schulten 1978): τ_c ~ ℏ / (Σ a²)^(1/2). This is purely algebraic — no MD needed if we use tabulated hyperfine constants for FAD.

**차박사:** I'll add the biophysical constraint: this is only meaningful for proteins where the radical pair separation distance is 10–25 Å (the Haberkorn kinetics window, Haberkorn 1976). Closer than 10 Å → spin exchange dominates; farther than 25 Å → singlet-triplet interconversion too slow. AlphaFold gives us the FAD binding site coordinates if the protein has an annotated flavin cofactor. We can check the distance between the two radical centers (FAD N5 and the Trp/Tyr terminal electron donor).

**노박사:** From biochemistry: let's be careful about overclaiming. The cryptochrome/photolyase family is well-defined (PF00627 + PF03441 Pfam domains). Outside this family, radical pair spin dynamics in proteins is speculative. My recommendation: **strict trigger — only activate this module when the protein has (a) FAD or FMN cofactor, AND (b) HMMER hit to PF00627 or PF03441 with E-value < 1e-10.** Don't apply it to generic flavoproteins.

**윤박사:** Agreed. The module name should be "Cryptochrome Radical Pair Analyzer" not "Radical Pair Spin-Dynamics Analyzer" — more honest about scope.

---

### 2C — Beratan-Onuchic Electron Tunneling Pathway (ETP) Mapper

**윤박사:** This is the most computationally tractable quantum biology feature. Beratan-Onuchic tunneling pathway theory (Beratan 1992, Science) computes the electron tunneling matrix element TDA through protein using a through-bond + through-space + through-hydrogen-bond model. Each pathway has a tunneling decay factor ε = Π(ε_bond) × Π(ε_space) × Π(ε_HB). No QM calculation required — it's a graph traversal problem on the protein structure.

**김박사:** This is implementable right now. Here's the algorithm: (1) Parse AlphaFold PDB → extract all heavy atoms. (2) Build a graph: nodes = atoms, edges = covalent bonds (distance < 1.8 Å) OR space jumps (1.8–4.0 Å) OR H-bonds (2.8–3.5 Å, donor-acceptor angle < 60°). (3) Assign decay factors: ε_bond = 0.6 per covalent bond, ε_space = 0.5 × exp(-β_space × r) with β_space = 1.7 Å⁻¹, ε_HB = 0.36. (4) Find top-k pathways from electron donor to acceptor using modified Dijkstra (maximize log(ε) = minimize -log(ε)). (5) Return pathway sequence + total TDA. Python implementation: `networkx` (already installed) + BioPython for PDB parsing.

**용박사:** I've implemented this before. The critical issue is defining donor and acceptor. For heme proteins: donor = Cys/His coordinating Fe, acceptor = surface Trp/Tyr. For Fe-S clusters: donor = Cys ligands, acceptor = next Fe-S cluster in the electron transfer chain. For flavoproteins: FAD N5 as acceptor, surface aromatic as donor. We need a **cofactor-based donor/acceptor assignment table** — not complicated, but requires careful curation. Also: only trigger this for redox-active proteins (heme, Fe-S, flavin, copper, molybdenum cofactors). Triggering on an antibody structure would produce meaningless pathways.

**차박사:** One implementation concern: AlphaFold coordinates have positional uncertainty (pLDDT). For ETP calculation, bond angles and distances must be accurate. My recommendation: only compute ETP for regions where pLDDT > 70. Flag residues in low-pLDDT regions as "pathway uncertainty high" if the computed tunneling path passes through them.

**Dr. Konberg:** From a drug design perspective, ETP pathways are fascinating because blocking a tunneling pathway (by binding at a key bridging residue) could selectively inhibit electron transfer in redox enzymes — a completely orthogonal mechanism to active site inhibition. The fragment hotspot map (already implemented) should cross-reference with ETP pathway residues. Residues that are both druggable hotspots AND ETP pathway nodes = very high-priority drug targets.

**Claudin:** This is the most immediately implementable quantum feature. Decision: **implement ETP mapper (networkx graph traversal, Beratan-Onuchic decay factors) as `analyzers/etp_mapper.py`. Trigger: redox cofactor present. Cross-reference with fragment_hotspot_map.**

---

### 2D — Vibrational Promoting Mode Profiler

**차박사:** The Schwartz-Schramm promoting vibration hypothesis (Schwartz 2009, Nature Chem Biol): certain low-frequency protein conformational modes (10–100 cm⁻¹) couple to the chemical step, reducing the effective barrier. These are "promoting vibrations" — correlated motions between remote residues and the active site that compress the donor-acceptor distance at the moment of transfer. We already have ANM/GNM conformational analysis implemented (ProDy). The question: can we identify promoting vibration candidates from ANM modes?

**용박사:** This is where I can contribute from computational chemistry. The Franck-Condon approach: a promoting vibration is one where the ANM mode displacement vector has significant projection onto the donor-acceptor distance vector. Algorithm: (1) Compute ANM modes (already done via `conformational_analyzer.py`). (2) For each mode k, compute the overlap: `Δr_DA · ê_k` where Δr_DA is the unit vector from donor to acceptor and ê_k is the mode displacement at donor/acceptor residues. (3) Modes with |overlap| > 0.3 are promoting vibration candidates. (4) Sort by overlap magnitude × frequency (lower frequency = softer = thermally accessible). This requires no new packages — just ProDy + the ETP mapper result to define donor/acceptor.

**김박사:** This is elegant because it chains two existing features: ETP mapper (identifies donor/acceptor) → vibrational profiler (uses ANM output + donor/acceptor) → promoting vibration candidates. Zero new dependencies. The output: list of (mode_number, frequency_cm-1, overlap_magnitude, residues_involved) with interpretation.

**윤박사:** The key citation here is Hammes-Schiffer & Benkovic 2006 (Science) — "Relating Enzyme Motion and Catalysis." And the Klinman lab's seminal work on coupled protein motions (Klinman & Kohen 2013, Annu Rev Biochem). We must cite these correctly in the code.

**Decision: implement as sub-module of `conformational_analyzer.py` — a new `identify_promoting_vibrations(anm_modes, donor_residue, acceptor_residue)` function. Trigger: only if ETP mapper result AND conformational ANM result both available.**

---

## Session 3 — ML Architecture Track

**Chair: 김박사**

**김박사:** Two items: contrastive PPI specificity predictor and LoRA family-specialized models. These are the most computationally intensive features proposed. Let me start with PPI.

---

### 3A — Contrastive PPI Specificity Predictor

**김박사:** Current PPI data comes from STRING (co-expression + experimental + text-mining aggregated score). This tells us *that* proteins interact — not *why* they interact specifically rather than with paralogs. The contrastive question: "EGFR interacts with GRB2 but not with GRB7 — what makes GRB2 binding specific?" This requires a model that encodes the binding interface, not just the sequence.

**차박사:** From biophysics: interface specificity is determined by shape complementarity (Sc score), electrostatic complementarity, and hot spot residues (Bogan & Thorn 1998, J Mol Biol). For known PPI structures (from PDB), these are computable. For predicted structures (AlphaFold-Multimer), they're less reliable but usable above pLDDT > 70 at the interface.

**용박사:** The contrastive learning approach: Siamese network where positive pairs are (protein A, its known interactor B) and negative pairs are (protein A, paralogs of B that don't interact). Feature vector: ESM-2 embedding of the binding interface residues (from STRING experimental edges + UniProt binding site annotation). Loss function: NT-Xent (normalized temperature-scaled cross entropy). Training data: all human PPI with experimental confirmation from STRING (score_type="experiments" > 700) filtered to pairs where both proteins have AlphaFold structures.

**김박사:** The problem is the training data pipeline. We'd need to: (1) fetch STRING experimental interactions for all human proteins, (2) extract ESM-2 embeddings for interface residues, (3) construct negative pairs (non-interacting paralogs), (4) train Siamese network. This is weeks of compute, not a single-session implementation. My recommendation: **defer full training to next month, but implement the inference pipeline now using a pre-trained PPI embedding model from HuggingFace.** There's `facebook/esm2_t33_650M_UR50D` for embeddings + a published contrastive PPI model: `Rostlab/prot_bert` variants fine-tuned on PPI. We use their checkpoint.

**노박사:** From biochemistry I want to add: the most important PPI specificity question for drug discovery is whether a small molecule can **selectively disrupt one interaction without affecting paralogs**. EGFR-GRB2 vs EGFR-SHC1 disruption selectivity. This is the "PPI selectivity landscape" — analogous to the kinase selectivity landscape we already have in `selectivity_analyzer.py`. We should design the PPI contrastive predictor output to feed directly into a "PPI disruption selectivity report."

**Decision: implement using pre-trained `facebook/esm2_t33_650M_UR50D` + cosine similarity at interface embeddings. Full contrastive training deferred to dedicated ML sprint. New file: `analyzers/ppi_specificity_analyzer.py`.**

---

### 3B — Family-Specialized LoRA Models

**김박사:** The proposal: fine-tune ESM-2 on specific protein families (kinases, GPCRs, ion channels, antibodies) using LoRA (Low-Rank Adaptation, Hu 2021). Family-specific models would outperform the general ESM-2 model on family-specific variant effect prediction.

**용박사:** LoRA is the right approach. For a 650M-parameter ESM-2, LoRA with rank r=16 adds only ~6M parameters per family. Fine-tuning on a V100 GPU (16GB VRAM) for ~2000 family sequences takes approximately 4 hours per family. Training data: (1) kinases — KinBase (Manning 2002) + KLIFS database, ~550 human kinase sequences + DMS data from MaveDB; (2) GPCRs — GPCRdb, ~800 sequences; (3) antibodies — OAS (Observed Antibody Space), millions of sequences.

**차박사:** The inference-time benefit: a LoRA-fine-tuned kinase model predicts drug resistance mutations 15-20% more accurately than base ESM-2 (demonstrated in Fraternali 2023 for kinase domain). For EGFR T790M resistance, the LoRA kinase model should correctly predict that T790M is tolerated by the kinase but reduces erlotinib binding — while base ESM-2 misclassifies it as structurally disruptive.

**Claudin:** From an implementation standpoint: the LoRA adapters are small (~50MB per family). We can host them on HuggingFace Hub under a `ProteinScope/` organization. At inference time, we detect protein family (from Pfam domain), download the relevant adapter, merge with base ESM-2. The `peft` library handles LoRA loading in 3 lines of code. The missing piece is training data curation and fine-tuning compute — that's a multi-week project.

**김박사:** Recommended families to train first, in order of therapeutic relevance:
1. **Kinases** — highest DMS dataset coverage (MaveDB), highest drug discovery relevance
2. **GPCRs** — 35% of all FDA-approved drugs target GPCRs
3. **Antibody VH/VL domains** — most important for biologic engineering
4. **Ion channels** — cardiac safety (hERG), CNS targets

**Decision: LoRA fine-tuning deferred — requires dedicated compute allocation. Immediate action: add `peft` to requirements.txt and implement the adapter-loading inference wrapper. When adapters are trained (external sprint), zero code changes needed — just upload adapters to HuggingFace Hub.**

---

## Session 4 — Engineering & UX Track

**Chair: Dr. Konberg**

**Dr. Konberg:** Four items: conformational MSM, APBS electrostatics, adaptive step skipping, and the two UX bugs. Let's move quickly.

---

### 4A — HMM/Markov State Models (Conformational Classification)

**차박사:** MSMs require MD trajectory data — typically 10–100 microseconds of simulation per protein. Without user-supplied trajectories, we cannot build meaningful MSMs. However: we can build **synthetic MSMs from NMA (Normal Mode Analysis) modes** using the approach of Berman & Bhargava (2021). ANM modes define a basin of harmonic motion; each mode defines a metastable state. The transition matrix between states can be estimated from mode frequencies (higher frequency → faster interconversion → higher transition rate). This is an approximation, but it's honest and zero-MD.

**용박사:** I've used this approach. The key limitation: NMA-based synthetic MSMs capture local harmonic basins but miss large-scale conformational switches (e.g., kinase DFG-in/DFG-out transition). For proteins with known conformational switches in the literature, we should supplement with PDB ensemble analysis — download all PDB structures for the target, cluster by RMSD, use cluster membership as pseudo-states.

**김박사:** PDB ensemble clustering is already feasible: BioPython's PDBParser + sklearn K-means on RMSD matrix (using `mdtraj` or just backbone distance matrix). If the protein has ≥ 5 PDB structures, we cluster them into 2–4 conformational states. If < 5 structures, we fall back to NMA synthetic MSM.

**Decision: implement PDB ensemble clustering (primary) + NMA synthetic MSM (fallback). New file: `analyzers/conformational_msm.py`. Trigger: always, with PDB structures fetched from RCSB PDB REST API.**

---

### 4B — APBS/DelPhi Electrostatics Engine

**차박사:** The APBS (Adaptive Poisson-Boltzmann Solver) computes the electrostatic surface potential of a protein — critical for understanding how charged drugs approach binding sites, for salt bridge analysis, and for proton binding (pKa). The APBS binary is large (~200MB) and requires PDB2PQR as a preprocessing step (assigns atom charges + radii).

**용박사:** Here's the feasibility path without requiring the user to install APBS locally: use the **APBS web service** (`https://server.poissonboltzmann.org/apbs`). It accepts PDB2PQR-processed files and returns the electrostatic potential in DX format. We then parse the DX grid and compute: surface potential average at binding site residues, salt bridge energies (Coulomb term), and isoelectric point agreement (compare computed vs UniProt pI).

**용박사:** However, the web service has 15-minute turnaround time — too slow for real-time SSE streaming. My alternative: use `pqr2pdb` + a lightweight Python Poisson-Boltzmann solver (`pbrd` package, CPU-only, 30-second runtime for proteins < 500 aa). Less accurate than APBS but sufficient for qualitative surface potential maps.

**노박사:** From biochemistry: the most clinically useful output from electrostatics isn't the full 3D potential map — it's **three numbers**: (1) active site electrostatic field strength (mV/Å) — relevant for charge-assisted catalysis; (2) binding pocket net charge at pH 7.4 — relevant for drug charge optimization; (3) most electrostatically unfavorable surface patch — relevant for antibody epitope engineering (de-immunization). Focus the implementation on these three, not the full potential map.

**Dr. Konberg:** Agreed. From medicinal chemistry: pocket charge at pH 7.4 determines whether a cationic drug (protonated amine, pKa ~9) binds electrostatically or sterically. This is why imatinib (piperazine pKa 7.7) has pH-dependent bioavailability — 30% higher plasma binding at pH 7.6 vs 7.2. Having this data in the protein report is genuinely novel for a platform like this.

**Decision: implement using PROPKA (already integrated in P4) for pKa + simplified Coulomb summation for pocket charge. Defer full APBS/pbrd to next sprint. New function: `analyzers/metalloenzyme_analyzer.py::compute_pocket_electrostatics()` — reuse existing PROPKA output.**

---

### 4C — Adaptive Step Skipping (Query-Aware Module Selection)

**김박사:** The core problem: currently `fetch_all()` runs all 40+ modules on every query. For "What is the 3D structure of TLR2?" we don't need PROTAC feasibility, AAV design, or CRISPR planning. This wastes 60 seconds and 40+ API calls on irrelevant analysis.

**Claudin:** My design proposal: add a `classify_query_modules()` call at the start of `fetch_all()` — Gemini classifies the query intent and returns a `modules: set[str]` of required module names. The `fetch_all()` function has a module registry mapping names to their trigger conditions and async function calls. Only the modules in `required_modules` are executed.

**용박사:** I want to add computational complexity classification. Some queries have compound intent: "Find PROTAC targets in the kinase domain of EGFR that have selectivity over HER2." This requires: druggability + PROTAC + selectivity + ETP (if EGFR has redox cofactor) + binding_energy. The module set is a *superset* of any single intent category. Gemini's classification must handle compound queries by returning a union of module sets.

**노박사:** The biochemist concern: module dependencies. If you skip the druggability module but include the PROTAC module, PROTAC will fail because it reads `druggability_report` from `ProteinRecord`. We need a **dependency graph** for modules. When a module is selected, all its dependencies are automatically included. Like a package manager: `npm install protac` also installs `druggability, binding_sites`.

**김박사:** Standard topological sort. The dependency graph is: `protac → druggability → fpocket → af_structure`. `selectivity → druggability → string_interactions`. `epistasis → esm2_scoring → sequence`. Easy to implement with a simple dict + DFS traversal.

**Dr. Konberg:** One more consideration: UI override toggles. Some users (bioinformatics experts) want to run specific modules manually. The UI should have a "Module Settings" panel where users can force-include or force-exclude modules regardless of query classification. This is a power-user feature for expert users — hide behind an "Advanced" toggle.

**Decision: implement query-intent classifier in `core/gemini_interpreter.py::classify_query_modules()` → module dependency graph in `core/module_registry.py` → DFS dependency resolution → pass `active_modules: set[str]` to `fetch_all()`. UI module toggles deferred to UX sprint.**

---

### 4D — UX Bugs (Items #9 and #10)

**Claudin:** Two UX issues from pre-reviewer feedback.

**Issue #9 — "Analyzing N proteins" banner has no click behavior.**

**Claudin:** Current: the banner "Analyzing 3 proteins: TLR2, IL1B, MMP9" is static text. Pre-reviewer asked: should clicking it expand a detail dropdown? My proposal: yes — clicking the banner expands a collapsible `<details>` panel showing per-protein status:
```
▾ Analyzing 3 proteins: TLR2, IL1B, MMP9
  ├── TLR2 ✓ completed in 12.3s — 14 databases, 8 analyzers
  ├── IL1B ⟳ running — step 3/6: Running IDP analysis...
  └── MMP9 ⏸ queued
```
After all proteins complete, the banner collapses to "✓ Analysis complete — 3 proteins analyzed."

**김박사:** Implementation: SSE events should carry per-protein metadata: `{"type": "step", "protein": "TLR2", "step": 3, "total_steps": 6, "elapsed_ms": 1234}`. The frontend maintains a `proteinStatus` dict and re-renders the banner on each SSE event. No backend changes to `fetch_all()` needed — just richer SSE payload.

**윤박사:** The quantum biologist input: in multi-protein comparison mode (e.g., "Compare TLR2 and TLR4"), the banner could show side-by-side cards. Not for this sprint — note for later.

**Issue #10 — Query lost on conversation turn; model restarts without restoring previous results.**

**Claudin:** Root cause hypothesis: on a new query submission, the SSE `EventSource` is recreated. The `loadReport()` function for the *previous* job is never called because the job status check only happens when the SSE closes with `type: "done"`. If the user types a new query before the previous job's `"done"` event fires, the previous job's report is lost from DOM.

**차박사:** Three-part fix: (1) Before starting a new SSE connection, check `localStorage` for a `last_job_id`. If it exists and is not the current job, call `loadReport(last_job_id)` to restore previous results (as collapsed cards). (2) Preserve the user's query text in `localStorage.setItem('last_query', queryText)` on form submit. On page load, if `last_query` exists, populate the input field. (3) When a new query is submitted, push a "New query started" divider card (timestamp + query text) to the chat, so the conversation history is visually clear.

**Dr. Konberg:** For the medicinal chemist use case: I often submit a query, read part of the result, then submit a refined query. If the first result disappears, I lose context. The collapsed cards approach (차박사's suggestion) is correct — previous results should persist as collapsed accordions, not disappear. Click to re-expand.

**Decision:**
- #9: Richer SSE payload with per-protein metadata + collapsible banner panel. Frontend change only.
- #10: `localStorage` preservation of `last_job_id` + `last_query` + restore-previous-results as collapsed cards. Frontend change only. Backend adds `"protein": gene_name` field to all `type: "step"` SSE events.

---

## Session 5 — Cross-Cutting Architectural Decisions

**Chair: Claudin**

### 5A — Citation Discipline

**노박사:** I want to formally propose a citation standard for this project. Every threshold, formula, and computed score must have a `# Citation: Author YEAR doi` comment. This already exists in some files but not consistently. Given this platform may be used in academic contexts, this is non-negotiable.

**Resolution: adopted as project invariant. Code reviewer (Dr. Konberg) will audit citation coverage before each release.**

---

### 5B — Confidence Propagation

**용박사:** Currently `EvidenceGrade` (experimental/computational/literature/ai_generated) is assigned per-module. But when a downstream module *reads* the output of an upstream module, it inherits the upstream uncertainty. Example: ETP mapper uses AlphaFold PDB (computational). Its output should be graded `computational` or lower — but if the tunneling prediction is then fed to Gemini synthesis, the output is `ai_generated` *based on* `computational` data. The compound uncertainty is worse than either individually.

**윤박사:** In quantum measurement theory, uncertainty propagates multiplicatively for independent sources and additively for correlated sources. For our system: we can define a simple compound grade:
```
experimental + experimental → experimental
experimental + computational → computational
computational + computational → computational (but note: double-inference)
computational + ai_generated → ai_generated
any + ai_generated → ai_generated
```
This is a partial order on evidence grades. The `DataProvenance` model should store `upstream_sources: List[DataProvenance]` to make the chain explicit.

**Resolution: add `upstream_sources: List[DataProvenance] = Field(default_factory=list)` to `DataProvenance` model. All derived modules populate this field.**

---

### 5C — Performance Budget

**김박사:** With 40+ modules and now potentially 50+ after this session, the median query time will exceed 3 minutes for a complex protein. We need a performance budget. My proposal:
- **Tier 1 (always run, < 30s)**: UniProt, sequence, AlphaFold URL, KEGG, STRING, ClinVar, drug interactions, Gemini summary
- **Tier 2 (run if triggered, < 60s each)**: ESM-2, IDP, PTM, metalloenzyme, conformational, host compatibility
- **Tier 3 (on-demand only, explicit user request)**: inverse folding, CRISPR, AAV, ETP mapper, LoRA models, APBS

**Claudin:** This aligns perfectly with the adaptive step skipping proposal. Tier 1 = always active; Tier 2 = activated by query intent; Tier 3 = explicitly requested by user. The module registry from Session 4C can encode the tier for each module.

**Resolution: module registry includes `tier: Literal[1, 2, 3]`. Query intent classifier activates Tier 2 modules as appropriate. Tier 3 requires explicit user action (button click or keyword in query).**

---

## Final Summary — Implementation Roadmap

### Immediate (< 2 weeks)

| Feature | File | Owner expertise |
|---|---|---|
| Disease Association Layer (OpenTargets + DisGeNET + HPO + ClinGen) | `fetchers/opentargets.py`, `fetchers/disgenet.py`, `fetchers/clingen.py`, `analyzers/disease_association_analyzer.py` | 김박사 (data integration) |
| ETP Mapper (Beratan-Onuchic, networkx) | `analyzers/etp_mapper.py` | 용박사 + 김박사 |
| Vibrational Promoting Mode (ANM overlap) | extend `analyzers/conformational_analyzer.py` | 차박사 + 용박사 |
| Conformational MSM (PDB ensemble + NMA fallback) | `analyzers/conformational_msm.py` | 차박사 |
| Pocket Electrostatics (PROPKA-based) | extend `analyzers/metalloenzyme_analyzer.py` | 용박사 + Dr. Konberg |
| Module Registry + Tier system | `core/module_registry.py` | 김박사 + Claudin |
| Query intent classifier | extend `core/gemini_interpreter.py` | Claudin |
| Module dependency graph + DFS resolver | `core/module_registry.py` | 김박사 |
| SSE per-protein metadata (richer payload) | `web/app.py` | Claudin |
| Frontend: collapsible banner + localStorage restore | `web/templates/index.html` | Claudin |
| `DataProvenance.upstream_sources` field | `core/evidence.py` | all |

### Next Sprint (2–6 weeks)

| Feature | Blocker |
|---|---|
| Cryptochrome Radical Pair Analyzer | Need HMMER subprocess + tunneling enzyme lookup table curation |
| Quantum Tunneling Flag (Tier 1 + 2) | Need tabulated ν‡ values for top 20 enzyme families |
| PPI Contrastive Specificity Predictor (inference only) | Need HuggingFace PPI embedding model selection |
| LoRA Kinase Adapter (inference wrapper) | Need `peft` + HuggingFace Hub adapter upload |
| Full APBS electrostatics | Need `pbrd` package evaluation |

### Deferred (dedicated ML/MD sprint)

| Feature | Why deferred |
|---|---|
| LoRA fine-tuning (kinase/GPCR/antibody/ion channel) | Weeks of compute + training data curation |
| Full contrastive PPI Siamese network training | Training infrastructure + negative pair curation |
| Full MSM from MD trajectories | Requires user-supplied MD data |

---

## Action Items

| # | Action | Responsible expertise | Target |
|---|---|---|---|
| A1 | Implement `fetchers/opentargets.py` (GraphQL, Ensembl ID from UniProt xrefs) | 김박사 domain | Week 1 |
| A2 | Implement `fetchers/disgenet.py` (REST, 1 req/sec rate limit) | 김박사 domain | Week 1 |
| A3 | Implement `fetchers/clingen.py` + `fetchers/hpo.py` | 김박사 domain | Week 1 |
| A4 | Implement `analyzers/disease_association_analyzer.py` (MONDO dedup + Gemini synthesis) | 노박사 + Claudin domain | Week 1–2 |
| A5 | Implement `analyzers/etp_mapper.py` (Beratan-Onuchic, networkx graph traversal) | 용박사 + 김박사 domain | Week 2 |
| A6 | Add `identify_promoting_vibrations()` to conformational_analyzer.py | 차박사 + 용박사 domain | Week 2 |
| A7 | Implement `analyzers/conformational_msm.py` (PDB cluster + NMA fallback) | 차박사 domain | Week 2–3 |
| A8 | Add `compute_pocket_electrostatics()` to metalloenzyme_analyzer.py | 용박사 + Dr. Konberg domain | Week 2–3 |
| A9 | Implement `core/module_registry.py` (tier + dependency graph + DFS) | 김박사 + Claudin domain | Week 1 |
| A10 | Add `classify_query_modules()` to gemini_interpreter.py | Claudin domain | Week 2 |
| A11 | Richer SSE step events (per-protein metadata) + frontend banner collapsible | Claudin domain | Week 1 |
| A12 | Frontend localStorage session restore + collapsed previous-result cards | Claudin domain | Week 1–2 |
| A13 | Add `upstream_sources` to DataProvenance; update all analyzers | all domain | Week 2–3 |
| A14 | Add `peft` to requirements.txt + LoRA inference wrapper stub | Claudin domain | Week 1 |

---

## Closing Remarks

**윤박사:** The Beratan-Onuchic ETP mapper is scientifically the most exciting feature we've discussed today. When ProteinScope can tell a researcher "here are the three electron tunneling pathways in your metalloprotein, and two of them pass through druggable fragment hotspots" — that's something no other publicly available platform does. I'm proud of where this is going.

**노박사:** From a biochemist's perspective: the disease association layer is the feature that transforms this from a "protein lookup tool" into a "protein-disease research engine." That's the difference between academic curiosity and clinical utility. The MONDO deduplication and ClinGen validity classification will be what distinguishes ProteinScope from DisGeNET's own interface.

**차박사:** The conformational MSM from PDB ensembles is underrated. Knowing that EGFR exists in three conformational states (DFG-in active, DFG-out inactive type I, αC-helix-out inactive type II) — and which drug binds which state — is fundamental to resistance mechanism prediction. We already have the conformational ANM analysis; PDB ensemble clustering is a minor extension with major clinical impact.

**용박사:** 오늘 회의에서 가장 인상적인 것은 모든 제안이 기존 infrastructure (networkx, BioPython, ProDy, PROPKA)를 활용한다는 점입니다. 새 heavy dependency가 거의 없습니다. 이건 engineering excellence의 증거입니다.

**김박사:** The module registry with tier classification is the architectural decision I'm most excited about. It transforms ProteinScope from a monolithic pipeline into a composable, query-aware system. That's the right abstraction for a platform that now has 50+ modules.

**Dr. Konberg:** One final note from medicinal chemistry: every feature we've designed today has a clear line to a drug discovery decision. ETP mapper → redox enzyme inhibition; disease association layer → target validation; pocket electrostatics → charge optimization; conformational MSM → selectivity by state. This is not academic tool-building. This is a medicinal chemistry decision-support system. That's what should be in the README.

**Claudin:** Closing the meeting. 14 action items assigned. Next meeting: after Week 3 implementation review. Minutes will be committed to the repository.

---

*Meeting adjourned — 2026-03-20*
*Minutes recorded by: Claudin, PhD (Claude Division)*
*Total duration: ~4 hours (simulated)*
*Next meeting: 2026-03-27 (tentative) — Implementation Review*
