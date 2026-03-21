# ProteinScope 그랜드 컨소시엄 회의록 — V4
## Codex GPT-5.4 × 3 + Gemini 2.5 Pro × 2 + Claude Sonnet 4.6 (의장)
## 2026-03-20

---

| 참석자 | 모델 | 역할 | 소속 |
|---|---|---|---|
| 김박사 | Codex GPT-5.4 | 전산 생물정보학자 | Stanford · MIT · Harvard |
| 차박사 | Codex GPT-5.4 | 생물물리학자 | Oxford |
| 용박사 | Codex GPT-5.4 | 전산 생물물리화학자 | SNU · KAIST |
| 노박사 | Gemini 2.5 Pro | 생화학자 | POSTECH |
| 윤박사 | Gemini 2.5 Pro | 양자생물학자 | CALTECH |
| 의장 | Claude Sonnet 4.6 | AI 바이오테크 개발자 | — |

---

## 의제

**연구 방향 쿼리(ROQ) 처리 아키텍처 설계**

사용자가 "말라세지아 모낭염에 효과적인 연고를 개발하려고 하는데 타겟이 될 만한 단백질은 무엇인가?"와 같이 질환/병원체는 알고 있지만 타겟 단백질을 모르는 경우, 기존 시스템은 `TaxonAmbiguityError`를 발생시켜 Saccharomyces, Penicillium, Agaricus 등 관계없는 종 선택 카드를 보여줬다. 이를 근본적으로 해결하는 아키텍처를 설계한다.

---

## 의장 (Claude) — 문제 정의

기존 흐름의 문제:
1. 사용자가 "말라세지아 모낭염"을 입력하면 UniProt 검색 실패
2. `TaxonAmbiguityError` 발화 → "곰팡이" 분류군 내 종 선택 카드 표시
3. 표시된 종: Saccharomyces cerevisiae, Penicillium chrysogenum, Agaricus bisporus, Neurospora crassa → **연구에 전혀 무관한 종들**

이것은 UI 버그가 아니라 **쿼리 라우팅 아키텍처의 근본적 결함**이다. 사용자는 단백질을 묻는 게 아니라 "어떤 단백질을 봐야 하는가"를 묻고 있다. 완전히 다른 처리 경로가 필요하다.

---

## 이박사 (Codex #1) — 탐지 우선순위 결정

`TaxonAmbiguityError`보다 먼저 ROQ를 잡아야 한다. 현재 `_resolve_accession()` 내 처리 순서는:

1. 경로/호환성/RNAi/역유전학 체크
2. NL 쿼리 파싱 (`parse_nl_protein_query`)
3. TaxonAmbiguityError ← **이 앞에 ROQ 체크 삽입 필요**

탐지 조건: **트리거 키워드** (연고, 치료, 개발, 항진균, 약물...) + **질환 키워드** (여드름, 감염, 피부염, folliculitis...) + **특정 단백질 이름 없음**

세 조건을 모두 만족할 때만 `TargetDiscoveryRequest` 예외를 발화한다. 이렇게 하면 "EGFR 항암 치료"처럼 단백질 이름이 명시된 경우는 기존 경로를 타게 된다.

---

## 정박사 (Codex #2) — 데이터 파이프라인 설계

5단계 파이프라인을 제안한다:

**Step 1 — Gemini 파싱**
자유 텍스트에서 `{organism_scientific, organism_common, disease_name, modality_hint, context_summary}` 추출.
한국어 입력도 처리 가능해야 한다.

**Step 2 — 프로테옴 fetch**
UniProt REST API: `reviewed:true AND organism_name:"{organism}"` 쿼리.
Reviewed 결과가 5개 미만이면 unreviewed 포함으로 fallback.
최대 80개 단백질 반환.

**Step 3 — 개별 단백질 스코어링**
- `essentiality_weight`: GO term 기반 (필수 대사/복제 = 1.0, 막 연관 = 0.7, 기타 = 0.4)
- `selectivity_bonus`: 인간 서열 동일성 < 30% → 1.5배, 30-60% → 1.0, > 60% → 0.3배
- `known_target_bonus`: 기존 알려진 항진균 타겟(ERG11, FKS1 등) + 0.3
- 최종 점수: `(essentiality × selectivity × literature) + known_bonus`

**Step 4 — 인간 호몰로그 검사**
`analyzers/sequence_aligner.py`의 기존 `align_pairwise()` 재사용.
인간 오솔로그를 UniProt에서 fetch → BLOSUM62 pairwise alignment → 동일성 % 계산.
동시 요청 Semaphore(3)으로 제한.

**Step 5 — Gemini 합성**
상위 10개 후보에 대해 기계적 근거 + 전체 치료 전략 생성.

---

## 이박사 (Codex #3) — 알려진 항진균 타겟 사전 구축

하드코딩된 타겟 사전이 필요하다. 스코어링 시 공식 문헌 기반 보너스를 부여한다:

```
Malassezia: LIP1, LIP2, LIP3, MGL1, FAS1, ERG11, ERG6, CYP51
Candida:    ERG11, FKS1, CDR1, MDR1, SAP1, SAP2, HWP1, ALS3
Aspergillus: CYP51A, FKS1, GEL1, CHS1, PksP, Lac1
Trichophyton: ERG11, SQS1, HMG1, CHS1, CHS2
Cryptococcus: ERG11, FKS1, LAC1, PKA1, CNA1
```

각 타겟에 대해 알려진 억제제 계열도 저장:
- ERG11 / CYP51 → 아졸계 (플루코나졸, 보리코나졸)
- FKS1 → 에키노칸딘계 (카스포펀긴, 미카펀긴)
- LIP1/2/3 → 지질분해효소 억제제 (올리스타트 유사체)

**인용**: Odds (2010) Antimicrob Agents Chemother; Rex (2016) Clin Infect Dis

---

## 노박사 (Gemini #1) — 데이터 모델 설계

`TargetDiscoveryReport`는 `ProteinRecord` 및 `models.py`에 포함시키지 않는다. 기존 `AntigenDiscoveryReport`와 동일한 독립 SSE 플로우로 처리.

```python
class TargetCandidate(BaseModel):
    rank: int
    accession: str
    gene_name: str
    protein_name: str
    organism: str
    target_score: float           # 0.0–1.0
    essentiality_note: str
    human_identity_pct: Optional[float]
    selectivity_risk: bool        # True if identity > 60%
    known_target: bool
    known_inhibitor_class: Optional[str]
    modality_hint: str
    gemini_rationale: str

class TargetDiscoveryReport(BaseModel):
    query: str
    disease_name: str
    organism_scientific: str
    organism_common: str
    candidates: List[TargetCandidate]
    context_summary: str
    gemini_strategy: str
    timestamp: datetime
```

---

## 윤박사 (Gemini #2) — 프론트엔드 UX 요구사항

연구자가 이 결과를 실제로 사용하려면:

1. **랭킹 테이블** — 점수 막대, 선택성 배지(초록/주황), 알려진 억제제 표시
2. **"Analyze →" 버튼** — 해당 단백질에 대한 전체 ProteinScope 분석을 즉시 실행
3. **맥락 요약** — 질환 배경 2문단 (접고 펼치기)
4. **치료 전략** — Gemini가 합성한 전체 전략 (접고 펼치기)
5. **면책 문구** — "계산 추정치이며 실험적 검증 필요"

---

## 합의 사항

| 항목 | 결정 |
|---|---|
| 탐지 시점 | `TaxonAmbiguityError` 이전, NL 파싱 이후 |
| 탐지 조건 | 트리거 키워드 + 질환 키워드 + 특정 유전자명 없음 (3중 조건) |
| 데이터 모델 | `models.py` 미포함, 독립 SSE 플로우 |
| 스코어 공식 | `(essentiality × selectivity × literature_proxy) + known_bonus` |
| 인간 호몰로그 | `align_pairwise()` 재사용, Semaphore(3) |
| 동시성 | asyncio.Semaphore(3) for UniProt fetch |
| 타겟 사전 | 5개 속(Genus) × 최대 8개 타겟 하드코딩 |
| 프론트엔드 | 랭킹 테이블 + "Analyze →" 버튼 + 접이식 전략 패널 |

---

## 검증 시나리오

| 쿼리 | 예상 결과 |
|---|---|
| `"말라세지아 모낭염 연고 개발하려는데 타겟 단백질은?"` | ERG11 (0.91), LIP1 (0.84), FAS1 (0.79) 상위 랭크 |
| `"Candida antifungal drug targets"` | ERG11, FKS1 최상위; CDR1은 efflux pump 주의 표시 |
| `"KRAS cancer therapy"` | 특정 유전자명 "KRAS" 존재 → ROQ 미발화 → 기존 단백질 분석 경로 |
| `"Malassezia globosa skin disease treatment"` | 영어 쿼리도 동일하게 처리 |

---

*회의록 작성: Claude Sonnet 4.6 · 2026-03-20*
