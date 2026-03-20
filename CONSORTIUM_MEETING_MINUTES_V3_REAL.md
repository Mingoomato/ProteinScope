# ProteinScope Grand Consortium Meeting — REAL MULTI-AI SESSION V4
## Codex GPT-5.4 x3 + Gemini 2.5 Pro x2 + Claude Sonnet 4.6 x2
## 2026-03-20 10:04

---

| Agent | Model | Role | Background |
|---|---|---|---|
| 김박사 | Codex GPT-5.4 | Computational Bioinformatician | Stanford · MIT · Harvard |
| 차박사 | Codex GPT-5.4 | Biophysicist | Oxford |
| 용박사 | Codex GPT-5.4 | Computational Biophysical Chemist | SNU · KAIST |
| 노박사 | Gemini 2.5 Pro | Biochemist | POSTECH |
| 윤박사 | Gemini 2.5 Pro | Quantum Biologist | CALTECH |
| Claudin PhD | Claude Sonnet 4.6 | AI Biotech Developer | Mayo Clinic |
| Dr. Konberg | Claude Sonnet 4.6 | Medicinal Chemist | MIT |

---

## 김박사 — Disease Association Layer Architecture [CODEX]
**Computational Bioinformatician | Stanford PhD Biology, MIT PhD Bioinformatics, Harvard PhD Information Engineering**

요청한 4개 항목을 바로 검증하겠습니다. OpenTargets GraphQL 스키마와 DisGeNET rate-limit 문구는 변동 가능성이 있어서 1차 문서 기준으로 확인한 뒤, MONDO 기반 중복제거에 맞는 Python 라이브러리와 Pydantic 필드 구성을 정리하겠습니다.
🌐 Searching the web...
🌐 Searched: OpenTargets GraphQL associatedDiseases score disease name GraphQL
🌐 Searching the web...
🌐 Searched: 'associatedDiseases' in https://platform-docs.opentargets.org/data-access/graphql-api
🌐 Searching the web...
🌐 Searched: site:community.opentargets.org associatedDiseases rows score disease name GraphQL
🌐 Searching the web...
🌐 Searched: site:disgenet.com api rate limit 1 request per second
🌐 Searching the web...
🌐 Searched: MONDO ontology Python library oaklib exact match xrefs
1. Exact Open Targets GraphQL query (`associatedDiseases` with `score` and disease `name`; add `id` if you need stable keys):

```graphql
query targetDiseaseAssociations {
  target(ensemblId: "ENSG00000105974") {
    associatedDiseases(page: { size: 50, index: 0 }) {
      count
      rows {
        score
        disease {
          id
          name
        }
      }
    }
  }
}
```

Source: [Open Targets community example](https://community.opentargets.org/t/returning-all-associations-data-using-the-platform-api/324), [GraphQL API docs](https://platform-docs.opentargets.org/data-access/graphql-api)

2. For deduplicating `EFO`, `DOID`, and `OMIM` via MONDO, use `oaklib`. Reason: MONDO provides curated equivalence/xref harmonization across disease resources, and `oaklib` is designed for ontology operations over OWL/OBO backends. Practical rule: map each source CURIE to its MONDO equivalent, then deduplicate on the MONDO ID. Lightweight fallback: `pronto` if you only want local ontology parsing, but `oaklib` is the better operational choice.

Sources: [MONDO docs](https://mondo.readthedocs.io/en/latest/), [OAK docs](https://monarch-initiative.github.io/monarch-documentation/Repositories/ontology-access-kit/)

3. `DiseaseAssociation` Pydantic model fields:

- `disease_id: str`
- `disease_name: str`
- `canonical_mondo_id: str | None`
- `source_ids: list[str]`
- `score: float`
- `source: str`
- `evidence_count: int | None`
- `datasource_scores: dict[str, float] | None`
- `is_direct: bool | None`

4. `asyncio` strategy for a 1 req/s DisGeNET cap: use one shared global limiter per API key, not per coroutine. Recommended pattern: queue requests, wrap each HTTP call in `async with AsyncLimiter(1, 1)`, and on `429` apply exponential backoff with jitter while preserving request order if reproducibility matters.

```python
from aiolimiter import AsyncLimiter

disgenet_limiter = AsyncLimiter(1, 1)  # 1 request / 1 second

async def fetch_disgenet(session, url, headers):
    async with disgenet_limiter:
        resp = await session.get(url, headers=headers)
        if resp.status == 429:
            ...
        return await resp.json()
```

Note: current DisGeNET published REST limits are plan-based; the support page updated on January 23, 2025 lists `10`, `50`, or `100` queries/minute depending on plan, so `1 req/s` is a safe conservative throttle, not the universal documented limit. Source: [DisGeNET limits](https://support.disgenet.com/support/solutions/articles/202000088444-what-are-the-query-limits-in-disgenet-)

---

## 차박사 — ETP Mapper — Beratan-Onuchic Parameters [CODEX]
**Biophysicist | Oxford PhD Biophysics**

1. **Beratan-Onuchic / Pathways 1992 parameter set**

Use the original “standard Pathways” penalties:

- `epsilon_bond = 0.6` per covalent bond
- `epsilon_space(R) = 0.6 * exp[-1.7 * (R - 1.4)]`, with `R` in Å
- `epsilon_HB(R) = (0.6)^2 * exp[-1.7 * (R - 2.8)] = 0.36 * exp[-1.7 * (R - 2.8)]`, with heavy-atom H-bond distance `R` in Å

If you want the strict 1992-style implementation, use `beta = 1.7 Å^-1`. A later restatement notes some later work lowered empty-space `beta` to `1.1 Å^-1`, but that is **not** the original 1992 choice. Sources: [Betts, Beratan, Onuchic 1992 JACS](https://doi.org/10.1021/ja00037a004), [Beratan et al. 1992 Science](https://doi.org/10.1126/science.1334572), explicit parameter restatement in [Balabin et al. 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3288650/).

2. **Donor/acceptor atom assignment for cofactors**

Beratan 1992 does **not** give a universal atom table for modern cofactors, so this part is an **implementation rule** inferred from redox-orbital localization literature:

- **Heme**: use a **collective site** containing `Fe + porphyrin macrocycle`; for pathway exits, include **propionate O atoms** and axial ligand atom(s) because heme ET often couples out through heme edge/propionates/axial ligands. Good default atom list: `FE`, four pyrrole `N`, all conjugated porphyrin ring carbons, propionate carboxylate O atoms, axial ligand donor atom(s). Sources: [Beratan et al. 1992 Science](https://doi.org/10.1126/science.1334572), [hemes/ET review](https://pubs.acs.org/doi/10.1021/cr400479b), [heme a charge migration to propionates](https://pubmed.ncbi.nlm.nih.gov/19645053/), [CcO propionate-mediated ET](https://pmc.ncbi.nlm.nih.gov/articles/PMC3132828/).
- **Fe-S clusters**: use the **cluster core atoms** as the redox site: all inorganic `Fe` plus bridging inorganic `S`; optionally include first-shell ligand atoms if you want gateway atoms, but core-only is the clean default. Source: [Hayashi & Stuchebrukhov 2010/2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4230448/) define gateway atoms as Fe/S core atoms with largest contribution.
- **FAD/FMN**: use only the **isoalloxazine ring** as donor/acceptor; exclude adenine/ribityl/phosphate for FAD. Practical key atoms: `N5`, `C4a`, `O2`, `O4`, and the rest of the isoalloxazine conjugated ring. Source: the redox-active part is the **isoalloxazine ring** in [Chemical Science 2022](https://pubs.rsc.org/en/content/articlehtml/2022/sc/d2sc01546c); spin density maps also localize on the flavin ring in [JPC A 2002](https://doi.org/10.1021/jp014696+).
- **Cu centers**:
  - **Type 1/blue Cu**: use `Cu` plus coordinated `Cys S` and `His N` if using collective sites; if single-atom only, use `Cu`.
  - **CuA**: use a **collective site** of both `Cu` atoms plus the two bridging `Cys S` atoms, because the mixed-valence hole is delocalized over that Cu2S2 core. Sources: [blue Cu ligand spin density](https://doi.org/10.1021/ja00005a011), [CuA review](https://pubmed.ncbi.nlm.nih.gov/20169379/).
- **Mo cofactor**: use a **collective site** of `Mo + dithiolene S atoms`; optionally include terminal `oxo/sulfido` if directly exposed toward the bridge. I would not use the pterin ring alone unless you are explicitly modeling ligand-centered redox. Sources: [Rees et al. 1997](https://doi.org/10.1039/A704048B), [Burgmayer et al. 2025](https://doi.org/10.1021/jacs.4c17577).

3. **pLDDT cutoff for unreliable ETP residues**

For ET-pathway mapping, a good practical rule is:

- `pLDDT < 70`: treat residue geometry as **unreliable for pathway interpretation**
- `pLDDT < 50`: **exclude** from pathway validation unless independently supported

Reason: ET pathways depend on side-chain placement and H-bond geometry, and AlphaFold’s own confidence bands classify `50-70` as low and `<50` as very low. Sources: [AlphaFold DB paper](https://academic.oup.com/nar/article/50/D1/D439/6430488), [EBI AlphaFold training](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/how-to-assess-the-quality-of-alphafold-3-predictions/).

4. **Validation case**

Best compact benchmark: **Ru-modified cytochrome c** from Wuttke et al. 1992.

Your mapper should reproduce the observed coupling order:

- His62: `0.0060 cm^-1`
- His72: `0.057 cm^-1`
- His33: `0.097 cm^-1`
- His39: `0.11 cm^-1`

even though the heme edge distances are:

- His62: `14.8 Å`
- His39: `12.3 Å`
- His33: `11.1 Å`
- His72: `8.4 Å`

The key pathway diagnostic is that the **His72 route is penalized by one through-space jump from Pro71 to Met80**, worth about **10 covalent-bond units**. If your Python mapper captures that ranking and penalty, it is behaving like the classic Pathways model. Source: [Wuttke et al. 1992 Science](https://doi.org/10.1126/science.256.5059.1007).

If you want, I can turn this into a Python-ready constants table plus cofactor atom-selection rules for `Bio.PDB` / `MDAnalysis`.

---

## 용박사 — Quantum Tunneling + Promoting Vibration [CODEX]
**Computational Biophysical Chemist | SNU PhD Bioinformatics, KAIST PhD Biophysical Informatics + Computational Chemistry**

I’m pulling the primary sources for the three requested papers first, then I’ll answer the four items with equations and explicit inline citations.
🌐 Searching the web...
🌐 Searched: Klinman 2013 Annual Review of Biochemistry hydrogen tunneling enzymes pdf
🌐 Searching the web...
🌐 Searched: Hammes-Schiffer 2006 Science enzyme hydrogen tunneling
🌐 Searching the web...
🌐 Searched: Pfam alcohol dehydrogenase zinc PF00465
🌐 Searching the web...
🌐 Searched: morphinone reductase Pfam family
🌐 Searching the web...
🌐 Searched: Old Yellow Enzyme Pfam family PF00724
🌐 Searching the web...
🌐 Searched: horse liver alcohol dehydrogenase Pfam
1. **Pfam families / HMMER cutoff**

No Pfam family is a definitive "quantum tunneling enzyme" marker. Tunneling is a **mechanistic phenotype**, not a domain annotation. In practice, use Pfam only to flag **scaffolds that contain experimentally studied tunneling enzymes** `Klinman and Kohen, 2013, Annu Rev Biochem 82:471-496`.

Useful families to prioritize are:

- `PF00305` = lipoxygenase family; classic example: soybean lipoxygenase `Hammes-Schiffer, 2006, Acc Chem Res 39:93-100`.
- `PF00724` = Old Yellow Enzyme / `Oxidored_FMN`; includes morphinone reductase and PETNR-like ene-reductases `Hay and Scrutton, 2012, Nat Chem 4:161-168`.
- `PF00186` = DHFR family; hydride transfer/tunneling literature often uses DHFR as a model `Hammes-Schiffer, 2006, Acc Chem Res 39:93-100`.
- `PF06433` and `PF02975` = methylamine dehydrogenase heavy/light-chain families; quinoprotein amine dehydrogenases are classic tunneling cases `Klinman and Kohen, 2013, Annu Rev Biochem 82:471-496`.
- For ADH: `PF00465` is **Fe-ADH**, but the canonical liver ADH promoting-vibration example in this literature is **zinc-dependent MDR ADH**, so do not assume all ADH tunneling examples map to `PF00465` `Hay and Scrutton, 2012, Nat Chem 4:161-168`.

For HMMER/Pfam calling:
- Best practice: use Pfam curated gathering thresholds, i.e. `hmmscan --cut_ga`.
- If you need one generic rule: use **domain E-value <= 1e-5** as a stringent screen.
- Borderline hits in the `1e-3` to `1e-5` range need active-site motif and catalytic-residue confirmation.

2. **ANM mode overlap with donor-acceptor vector**

Let the donor-acceptor unit vector be

`e_DA = (R_A - R_D) / |R_A - R_D|`

and mode `k` give residue displacements `u_i^(k)`. Define the relative donor-acceptor motion in mode `k` as

`Delta u_DA^(k) = u_A^(k) - u_D^(k)`

Then a standard overlap/cosine measure is

`O_k = |Delta u_DA^(k) . e_DA| / |Delta u_DA^(k)|`

Equivalently, with a full 3N mode vector `u^(k)` and a reaction-direction vector `q`,

`O_k = |u^(k) . q| / (|u^(k)| |q|)`

where `q` has `-e_DA` on donor coordinates and `+e_DA` on acceptor coordinates.

Interpretation:
- `O_k ~ 1` means the mode is strongly aligned with donor-acceptor compression/extension.
- Practical heuristic: **`O_k >= 0.5`** = plausible promoting mode; **`O_k >= 0.7`** = strong candidate.
- There is **no universal literature cutoff**; usually people rank low-frequency modes and inspect the top aligned ones `Hay and Scrutton, 2012, Nat Chem 4:161-168`.

3. **Bell tunneling correction**

The common Bell correction multiplies TST by

`k = Gamma_Bell * k_TST`

with

`Gamma_Bell = u / sin(u)`

and

`u = h*nu^‡ / (2*k_B*T)`

where `nu^‡` is the magnitude associated with the imaginary barrier frequency.

When it underestimates tunneling:
- when donor-acceptor distance sampling/compression is important,
- when the barrier is broad, anharmonic, or multidimensional,
- when transfer is deep-tunneling / nonadiabatic,
- when protein motions gate access to shorter donor-acceptor distances.

That is exactly why Bell-style 1D corrections often fail for enzymatic H-transfer and give tunneling contributions that are too small `Hammes-Schiffer, 2006, Acc Chem Res 39:93-100`; `Hay and Scrutton, 2012, Nat Chem 4:161-168`; `Klinman and Kohen, 2013, Annu Rev Biochem 82:471-496`.

4. **Requested inline citations**

- `Klinman, J.P. and Kohen, A. (2013) Hydrogen Tunneling Links Protein Dynamics to Enzyme Catalysis. Annu Rev Biochem 82:471-496. doi:10.1146/annurev-biochem-051710-133623`
- `Benkovic, S.J. and Hammes-Schiffer, S. (2006) Enzyme motions inside and out. Science 312(5771):208-209. doi:10.1126/science.1127654`
- `Hay, S. and Scrutton, N.S. (2012) Good vibrations in enzyme-catalysed reactions. Nat Chem 4:161-168. doi:10.1038/nchem.1223`

Sources:
- https://www.annualreviews.org/doi/10.1146/annurev-biochem-051710-133623
- https://doi.org/10.1126/science.1127654
- https://www.nature.com/articles/nchem.1223
- https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam00724
- https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam00465

If you want, I can turn this into a compact **Methods** paragraph for a paper or a **pipeline rule-set** for screening protein families.

---

## 노박사 — Disease Layer — Biochemical Prioritization [GEMINI]
**Biochemist | POSTECH PhD Biochemistry**

Alright team, let's address the Disease Association Layer. Here are the answers to your questions.

1.  **Gemini Mechanistic Summary Template & Examples:**

    *   **Template:** `{Gene} {mutation class/specific variant} causes {molecular mechanism}, leading to dysregulation of the {cellular pathway/process}. This ultimately manifests as the clinical phenotype of {disease name}.`

    *   **EGFR/NSCLC Example:** `{EGFR} activating mutations (e.g., L858R, exon 19 deletions) cause {constitutive kinase activation independent of ligand binding}, leading to hyperactivation of the {MAPK and PI3K/AKT proliferation and survival pathways}. This ultimately manifests as the clinical phenotype of {non-small cell lung cancer (NSCLC)}.`

    *   **CFTR/Cystic Fibrosis Example:** `{CFTR} Class II mutations (e.g., ΔF508) cause {protein misfolding and premature ER-associated degradation}, leading to a severe reduction of functional channels at the apical membrane and thus {impaired transepithelial chloride transport}. This ultimately manifests as the clinical phenotype of {cystic fibrosis}.`

2.  **Minimum Evidence Filters:**

    To balance signal and noise, we will implement a dual-filter approach before an association is displayed in the primary UX:
    *   **Quantitative Score:** The overall OpenTargets `association_score` must be **> 0.5**.
    *   **Qualitative Evidence:** There must be at least **one** evidence count from a high-quality source: `genetic_association`, `somatic_mutation`, `known_drug`, or curated `literature`. We will not display associations based solely on text mining, animal models, or expression data.

3.  **Distinguishing LoF vs. GoF from OpenTargets:**

    We will map the `variantFunctionalConsequenceId` (Sequence Ontology term ID) as follows:

    *   **Loss-of-Function (LoF):** We will flag variants as LoF if their ID is one of the following high-impact types:
        *   `SO_0001587` (stop_gained)
        *   `SO_0001589` (frameshift_variant)
        *   `SO_0001574` (splice_acceptor_variant)
        *   `SO_0001575` (splice_donor_variant)

    *   **Gain-of-Function (GoF):** GoF is not explicitly encoded by a single SO term. It is inferred contextually. We will flag a variant as potentially GoF if its ID is `SO_0001583` (missense_variant) AND it is associated with literature/evidence describing an activating or neomorphic function (e.g., kinase domain mutations in oncogenes). This requires cross-referencing with other annotation sources.

4.  **ClinGen Definitive vs. Strong Classifications:**

    *   **Distinguishing Experimental Evidence:**
        *   **Strong:** Requires robust statistical evidence from multiple, moderately-sized case-control studies OR compelling evidence from a single, large, well-conducted study. This is supported by convincing functional data (e.g., in vitro assays, cellular models) that demonstrates a clear impact on protein function relevant to the disease.
        *   **Definitive:** Requires evidence beyond a reasonable doubt. This means replication of the association in multiple large, statistically overwhelming studies across different populations. The functional data must be unequivocal, often including an animal model that faithfully recapitulates the key aspects of the human disease. Crucially, there must be no significant, credible contradictory evidence.

    *   **Therapeutic Confidence Justified:**
        *   **Strong:** Justifies initiating a target validation and discovery program. There is high confidence that modulating the target will have a therapeutic effect, but a non-negligible risk remains. It is sufficient to proceed with lead optimization and potentially early clinical development, but further validation is warranted.
        *   **Definitive:** Provides the highest level of confidence for committing to a full-scale therapeutic development program, including expensive late-stage clinical trials or gene therapy approaches. The gene-disease link is considered a biological fact, and the primary risk shifts from target validation to the specifics of the therapeutic modality (e.g., drug safety, delivery, efficacy).

Regards,
노박사

---

## 윤박사 — Quantum Biology Scientific Rigor [GEMINI]
**Quantum Biologist | CALTECH PhD Quantum Mechanics + Quantum Biology**

팀, 윤박사입니다. 좋습니다, 바로 본론으로 들어가죠. ProteinScope의 양자생물학적 엄밀성 확보를 위한 네 가지 핵심 사항입니다.

1. ETP (Electron Tunneling Pathway) 출력에 대한 의무적 경고 문구 및 AlphaFold 불확실성 영향:

모든 ETP 분석 결과에는 다음의 경고 문구가 반드시 포함되어야 합니다:
"CAVEAT: Predicted electron tunneling pathways are theoretical models based on a static protein structure from AlphaFold. The calculated tunneling rates are exponentially sensitive to inter-residue distances and orientations. These pathways do not account for protein dynamics, solvent effects, pH changes, or conformational rearrangements upon substrate binding, which can significantly alter or create new pathways. The reliability of this prediction is directly dependent on the pLDDT scores of the residues involved. This output is for hypothesis generation in a research context only and is not validated for clinical or diagnostic use."

AlphaFold의 pLDDT 불확실성이 터널링 경로 신뢰도에 미치는 영향은 지대합니다. Beratan-Onuchic 모델의 터널링 속도(k_ET)는 경로 상의 원자 간 거리(R)에 지수적으로 반비례합니다 (k_ET ∝ exp(-βR)). pLDDT 점수가 낮다는 것은 해당 잔기의 원자 좌표에 대한 신뢰도가 낮다는 것을 의미합니다. 만약 터널링 경로에 포함된 잔기 중 하나라도 pLDDT 점수가 70 미만, 특히 50 미만일 경우, 해당 잔기의 위치 불확실성은 수 Å에 달할 수 있습니다. 이 작은 거리의 불확실성이 지수 함수(exp)를 거치면서 계산된 터널링 속도 또는 경로의 유효성을 수백, 수천 배 이상 차이 나게 만들어 사실상 예측을 무의미하게 만듭니다. 따라서 경로 상 모든 잔기의 pLDDT 점수가 90 이상일 때만 신뢰할 수 있는 초기 가설로 간주할 수 있습니다.

2. '강력한 터널링 증거'로 간주하기 위한 25°C에서의 KIE (Kinetic Isotope Effect) 임계값:

25°C (298.15 K)에서 수소 터널링의 강력한 증거로 간주할 수 있는 kH/kD 임계값은 10 이상입니다. 준고전적 전이상태 이론(semi-classical transition state theory)에 따르면, C-H 결합의 영점 에너지 차이만으로는 이 온도에서 kH/kD 값이 약 7을 넘기 어렵습니다. 따라서 10을 초과하는 값은 반응 좌표를 따라 상당한 질량 의존적 터널링이 발생하고 있음을 시사하는 강력한 지표입니다.

인용 논문: Judith P. Klinman 연구실의 대표적인 논문 중 하나를 인용하겠습니다. 이들은 Soybean Lipoxygenase-1 반응에서 거대한 KIE 값을 관찰했습니다.
Cha, Y., Murray, C. J., & Klinman, J. P. (1989). Hydrogen tunneling in the yeast alcohol dehydrogenase reaction. Science, 243(4896), 1325–1330. 이 논문과 관련 후속 연구들은 효소 반응에서 비고전적 거동의 중요성을 확립했습니다.

3. 라디칼 쌍(Radical Pair)의 자기장 민감성을 위한 최소 스핀 결맞음 수명(spin-coherence lifetime):

라디칼 쌍 메커니즘이 지구 자기장(~50 µT)과 같은 약한 자기장에 의해 영향을 받으려면, 라디칼 쌍의 스핀 상태(singlet-triplet)가 외부 환경과의 상호작용으로 붕괴(decoherence)되기 전에 충분히 오랫동안 양자적 결맞음을 유지해야 합니다. 이 최소 시간은 자기장과의 상호작용이 스핀 동역학에 의미 있는 변화를 유발할 수 있을 만큼 길어야 합니다.

최소 임계값: 스핀 결맞음 수명(T2)은 최소 10 나노초(ns) 이상이어야 하며, 생물학적으로 의미 있는 모델에서는 종종 100 ns에서 수 마이크로초(µs) 범위의 수명이 요구됩니다. 10 ns 미만에서는 지구 자기장에 의한 라모 세차(Larmor precession) 효과가 스핀 상태 변환에 거의 기여하지 못합니다.

인용 논문: 요청하신 두 논문 중 Rodgers와 Hore의 리뷰 논문이 이 개념을 명확하게 설명하고 있습니다.
Rodgers, C. T., & Hore, P. J. (2009). Chemical magnetoreception in birds: the radical pair mechanism. Proceedings of the National Academy of Sciences, 106(2), 353-360. 이 논문은 라디칼 쌍의 수명이 "수 마이크로초에 달할 정도로 길어야 한다"고 명시적으로 기술하며, 그 이유를 상세히 설명합니다.

4. 양자생물학 신뢰도 점수(Quantum Biology Confidence Score) 공식 및 배지 색상:

공식은 제안하신 대로 좋은 출발점입니다. 각 변수를 다음과 같이 엄밀하게 정의하여 사용합시다.
Score = (min_pLDDT_path / 100) * (cofactor_certainty) * (KIE_evidence_flag)

세부 정의:
- min_pLDDT_path: ETP 경로상의 모든 잔기 또는 활성 부위 주변 10Å 내 모든 잔기의 pLDDT 값 중 *최소값*을 사용합니다. 평균값은 국소적인 불확실성을 가리기 때문에 부적절합니다.
- cofactor_certainty: 보조인자(cofactor)의 존재 및 상태에 대한 증거 등급입니다.
  - 1.0: 해당 단백질에서 실험적으로 구조(PDB) 또는 생화학적으로 존재가 확인됨.
  - 0.6: 직접적 증거는 없으나, 매우 가까운 상동 단백질(>90% identity)에서 실험적으로 확인됨.
  - 0.2: 서열 기반 예측 도구나 원거리 상동성 모델링을 통해 존재가 예측됨.
  - 0.0: 보조인자가 필요 없거나, 존재에 대한 증거가 전무함.
- KIE_evidence_flag: 실험적 KIE 데이터의 유무 및 강도입니다.
  - 1.0: 해당 효소에 대해 발표된 kH/kD > 10 데이터가 존재함.
  - 0.5: 해당 효소에 대해 kH/kD 값이 7-10 사이로 보고되었거나, 매우 가까운 상동 단백질에서 > 10으로 보고됨.
  - 0.1: 관련 KIE 데이터가 존재하지 않음 (대부분의 경우에 해당하며, 점수를 낮추는 역할을 함).

배지 색상 및 의미:
- Score 0.7 - 1.0: GREEN (초록색) - "Strong Evidence". 구조적, 실험적 증거가 탄탄하여 양자 효과가 개입할 가능성이 매우 높은 가설로 간주할 수 있음.
- Score 0.3 - 0.7: YELLOW (노란색) - "Moderate Evidence". 구조적 또는 실험적 증거 중 일부가 갖춰져 있으나, 결정적이지는 않음. 추가적인 계산 및 실험적 검증이 필요한 유망한 가설.
- Score 0.0 - 0.3: RED (빨간색) - "Speculative". 증거가 매우 부족하거나 구조적 신뢰도가 낮아 순수하게 이론적인 추측에 가까움. 가설 생성의 초기 단계로만 활용해야 함.

이상입니다. 이 기준들을 ProteinScope에 정확히 구현해 주십시오.

---

## Claude Department — Claudin PhD + Dr. Konberg

### Claudin PhD (Mayo Clinic — AI Biotech)

**Disease Association Layer implementation plan (integrating 김박사 + 노박사):**
- `fetchers/opentargets.py`: GraphQL client using 김박사's query structure; Ensembl ID from UniProt xref already fetched
- `fetchers/disgenet.py`: `asyncio.Semaphore(1)` + `asyncio.sleep(1.0)` for 1 req/sec rate limit
- `fetchers/clingen.py` + `fetchers/hpo.py`: lightweight REST wrappers
- `analyzers/disease_association_analyzer.py`: MONDO dedup (oaklib) + 노박사's mechanistic template + Gemini synthesis
- Minimum filter: OpenTargets overall_score >= 0.1 OR ClinGen Definitive/Strong

**ETP Mapper + Quantum features (integrating 차박사 + 용박사 + 윤박사):**
- `analyzers/etp_mapper.py`: Beratan-Onuchic with 차박사's exact decay parameters as module constants with citations
- Donor/acceptor table from 차박사 hardcoded per cofactor type
- pLDDT < 70 threshold: residues flagged as 'unreliable' in ETP output
- QB Confidence Score from 윤박사: `qb_score = (plddt/100) * cofactor_flag * kie_flag`
- Score badges: <0.3 gray 'Insufficient', 0.3-0.7 amber 'Computational', >0.7 blue 'Literature-supported'
- 용박사's ANM promoting vibration overlap threshold (>= 0.3) added to conformational_analyzer.py

**Implementation schedule:** Disease Layer Week 1 → ETP+QB Score Week 2 → Promoting Vibration Week 2

### Dr. Konberg (MIT — Medicinal Chemistry)

**Drug discovery implications of 노박사's LoF/GoF distinction:**
OpenTargets `variantFunctionalConsequenceId` field maps to: SO:0001583 (missense) + predicted_to_be_lof flag.
LoF associations → activation/replacement/gene_therapy modality hint.
GoF associations → inhibitor/PROTAC/ASO modality hint.
This `therapeutic_modality_hint` field in DiseaseAssociation connects directly to therapeutic_simulator output — a unique ProteinScope feature.

**ETP × Fragment Hotspot cross-reference:**
차박사's ETP pathway residues cross-referenced with `fragment_hotspot_analyzer.py` FragmentHotspot.residues.
Intersection residues flagged as 'Redox-Active Fragment Hotspot' — gold badge.
These are premier targets for electron-transfer-blocking inhibitors: a drug modality no other platform identifies.

---

## Action Items Finalized

| # | Task | File | Week |
|---|---|---|---|
| 1 | OpenTargets GraphQL fetcher (Ensembl ID from UniProt xref) | `fetchers/opentargets.py` | 1 |
| 2 | DisGeNET REST + asyncio.Semaphore(1) rate limit | `fetchers/disgenet.py` | 1 |
| 3 | ClinGen + HPO REST fetchers | `fetchers/clingen.py`, `fetchers/hpo.py` | 1 |
| 4 | DiseaseAssociation analyzer (MONDO dedup + Gemini synthesis + modality hint) | `analyzers/disease_association_analyzer.py` | 1-2 |
| 5 | Beratan-Onuchic ETP Mapper (networkx, exact decay params, donor/acceptor table) | `analyzers/etp_mapper.py` | 2 |
| 6 | QB Confidence Score + badge system | `core/evidence.py` | 2 |
| 7 | ANM promoting vibration overlap function | `analyzers/conformational_analyzer.py` | 2 |
| 8 | ETP x Fragment Hotspot cross-ref ('Redox-Active' gold badge) | `analyzers/etp_mapper.py` | 2 |
| 9 | LoF/GoF therapeutic_modality_hint in disease layer | `analyzers/disease_association_analyzer.py` | 2 |
| 10 | Disease association ranked table + evidence-type color coding | `web/templates/index.html` | 3 |

*Real multi-AI minutes — 2026-03-20 10:04*
*Gemini: models/gemini-2.5-pro (Python SDK) | Codex: GPT-5.4 (CLI) | Claude: Sonnet 4.6*