# ProteinScope 그랜드 컨소시엄 회의록 — V5
## Codex GPT-5.4 × 3 + Gemini 2.5 Pro × 2 + Claude Sonnet 4.6 (의장)
## 주제: 상용 시스템 벤치마킹 및 뷰어 API 통합 전략
## 2026-03-20

---

| 참석자 | 모델 | 역할 | 소속 |
|---|---|---|---|
| 이박사 | Codex GPT-5.4 | UI/UX 전문가 | — |
| 김박사 | Codex GPT-5.4 | 바이오인포매틱스 DB 전문가 | — |
| 정박사 | Codex GPT-5.4 | 프론트엔드 아키텍처 전문가 | — |
| 나박사 | Gemini 2.5 Pro | 연구자 UX 전문가 | — |
| 오박사 | Gemini 2.5 Pro | 백엔드 안정성 전문가 | — |
| 의장 | Claude Sonnet 4.6 | AI 바이오테크 개발자 | — |

---

## 의제

**사용자 피드백**: "실제 우리 시스템에서 진행하는 분석의 결과가 너무 허접해. 실제 상용화된 시스템의 결과 도출 화면을 찾아보고 벤치마킹을 해서 가져오거나, 실제로 상용화된 시스템이 온라인 뷰어 API를 지원하는 것이 있다면 우리 Tool에 삽입하는 것이 좋은 방향인 것 같아."

---

## 의장 (Claude) — 현황 진단

현재 ProteinScope 출력의 세 가지 핵심 문제:

1. **구조 시각화**: Mol.js 뷰어가 AlphaFold PDB만 로드하며 실험 결정 구조(X-ray/cryo-EM) 무시
2. **발현 데이터**: 조직별 발현이 단순 텍스트 나열 — ProteomicsDB, GTEx 히트맵과 비교하면 정보 밀도 1/10 수준
3. **약물-단백질 상호작용**: ChEMBL은 화학구조 이미지 + IC50 분포 + assay 분류를 동시 제공 — 우리 카드는 이름과 기전만 나열

---

## 이박사 (Codex #1) — 상용 시스템 현황 진단

현재 최고 수준의 단백질 분석 플랫폼을 기준으로 우리 시스템 격차를 분석했다:

**UniProt 웹사이트**
- 특징: 도메인/PTM/변이체/구조를 인터랙티브 트랙으로 표시
- 우리 현황: 텍스트 목록

**STRING DB**
- 특징: 증거 채널별 색상 구분 PPI 네트워크, 신뢰도 슬라이더
- 우리 현황: Cytoscape 그래프 (유사하나 시각적 품질 낮음)

**ChEMBL**
- 특징: 2D 화학구조 이미지 + IC50/Ki 분포 + assay 타입 분류 동시 표시
- 우리 현황: 약물명 + 기전 텍스트만

**AlphaFold DB**
- 특징: pLDDT 잔기별 바 차트 + 4색 신뢰도 색상 체계 표준화
- 우리 현황: 단일 평균 숫자

**ProteomicsDB**
- 특징: 조직 × 세포타입 발현 히트맵
- 우리 현황: 없음

---

## 김박사 (Codex #2) — API 지원 현황 조사

상용 시스템별 무료 API 지원 현황:

### 🟢 즉시 통합 가능 (무료, 안정적)

| 시스템 | 제공 내용 | API 방식 |
|---|---|---|
| **RCSB PDB** | X-ray/cryo-EM 결정 구조 3D 뷰어 | NGL.js 직접 로드 |
| **ChEMBL** | 화합물 2D 구조 SVG 이미지 | `https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}.svg` |
| **EBI ProTVista** | UniProt 공식 feature annotation 트랙 | Web Component (`protvista-uniprot`) |
| **AlphaFold DB** | 잔기별 pLDDT 신뢰도 배열 | `GET /api/prediction/{uniprot_id}` |
| **STRING DB** | PPI 네트워크 이미지 | `https://string-db.org/api/image/network?identifiers={gene}&species={taxon}` |

### 🟡 조건부 통합 가능

| 시스템 | 조건 |
|---|---|
| Human Protein Atlas | 무료, CC BY-SA 라이선스 표기 필요 |
| ProteomicsDB | 무료 REST API, 장기 안정성 불확실 |
| OMIM | 무료 API key 신청 필요 |

### 🔴 통합 비권장

| 시스템 | 이유 |
|---|---|
| Schrödinger Glide | 상용 라이선스 |
| PDB-REDO | 스키마 변경 잦음 |
| I-TASSER | 서버 제한 심각 |

---

## 정박사 (Codex #3) — 기술 구현 패턴 분류

세 가지 통합 패턴으로 분류:

**패턴 1 — 이미지 URL 임베드** (가장 쉬움, 즉시 구현 가능)
```html
<!-- ChEMBL 화합물 2D 구조 -->
<img src="https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL25.svg"
     style="width:120px; height:100px;">

<!-- STRING PPI 네트워크 이미지 -->
<img src="https://string-db.org/api/image/network?identifiers=EGFR&species=9606&limit=10">
```
CORS 문제 없음 (이미지 URL은 same-origin 정책 미적용).

**패턴 2 — Web Component 임베드** (중간 난이도, 최고 품질)
```html
<script type="module"
  src="https://cdn.jsdelivr.net/npm/protvista-uniprot/dist/protvista-uniprot.js">
</script>
<protvista-uniprot accession="P00533"></protvista-uniprot>
```
EBI 공식 컴포넌트. `accession` attribute만 전달하면 EBI 서버에서 모든 데이터를 자동 fetch하여 UniProt 웹사이트와 동일한 트랙을 렌더링함.

**패턴 3 — 클라이언트 사이드 API fetch + D3 렌더링**
```javascript
const resp = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);
const data = await resp.json();
// confidenceScore 배열 → D3 바 차트
```

---

## 나박사 (Gemini #1) — 연구자 관점 우선순위 결정

연구자가 새 도구를 처음 열었을 때 "이건 제대로 된 도구다"라고 느끼게 하는 요소 TOP 3:

**1위: ProTVista Feature Viewer (EBI Web Component)**
- 이유: 연구자들이 UniProt에서 매일 보는 UI와 동일 → 즉시 신뢰감 형성
- 구현 난이도: 한 줄 코드로 완성
- 한 줄 `<protvista-uniprot accession="...">` 태그가 UniProt 공식 사이트와 동일한 feature track 렌더링

**2위: ChEMBL 화합물 구조 이미지**
- 이유: 텍스트 약물 목록이 시각적 약물 카드로 즉시 전환
- 구현 난이도: `<img>` 태그 한 줄 추가
- ChEMBL ID는 이미 기존 ChEMBL fetcher 출력에 포함

**3위: AlphaFold pLDDT 신뢰도 바 차트**
- 이유: 단일 평균 숫자보다 잔기별 분포가 훨씬 많은 정보 전달
- EBI 공식 4색 체계: >90=진파랑, 70-90=시안, 50-70=노랑, <50=주황
- AlphaFold DB API에서 `confidenceScore` 배열 무료 제공

---

## 오박사 (Gemini #2) — 백엔드 안정성 검토

**외부 API 의존성 리스크 분석:**

ProTVista, STRING 이미지, ChEMBL 이미지는 모두 **클라이언트 사이드** 요청이다. 백엔드 코드 변경 없이 `index.html`에만 추가하면 된다. 해당 API가 다운되어도 해당 섹션만 빈칸 표시 — graceful degradation 자동 처리.

**STRING API 주의사항:**
`identifiers` 파라미터에 특수문자 포함 시 `encodeURIComponent()` 처리 필수.

**AlphaFold DB API 주의사항:**
응답이 배열 형태 `[{confidenceScore: [...]}]`. 장 길이가 500 이상이면 다운샘플링 필요.

**오류 처리 패턴:**
```html
<img src="..." onerror="this.style.display='none'">
```
이미지 로드 실패 시 요소 자체를 숨김 — 사용자에게 오류 노출 없음.

---

## 합의 사항 및 우선순위

| 기능 | 구현 방식 | 백엔드 변경 | 우선순위 |
|---|---|---|---|
| ProTVista feature annotation 트랙 | EBI Web Component (CDN) | 불필요 | **즉시 (P1)** |
| ChEMBL 2D 화합물 구조 이미지 | `<img src="...svg">` | 불필요 | **즉시 (P1)** |
| AlphaFold pLDDT 잔기별 바 차트 | AlphaFold DB API + D3 | 불필요 | **즉시 (P1)** |
| STRING PPI 네트워크 이미지 | `<img src="...">` | 불필요 | **즉시 (P1)** |
| Human Protein Atlas 조직 발현 히트맵 | `fetchers/protein_atlas.py` 확장 + D3 | 필요 | 1주 내 (P2) |
| RCSB 결정 구조 NGL 뷰어 | NGL.js 신규 로드 | 불필요 | 1주 내 (P2) |

**P1 기능의 공통 특징**: 모두 `index.html` 수정만으로 완성. 서버 재시작 불필요. 리스크 제로.

---

## 구현 세부 결정

### ProTVista
- CDN: `https://cdn.jsdelivr.net/npm/protvista-uniprot/dist/protvista-uniprot.js`
- `<head>`에 `type="module"` 스크립트 태그 추가
- `loadReport()` 내 `setTimeout(() => showProTVistaFeatureViewer(uniprot_id), 400)` 패턴

### AlphaFold pLDDT 차트
- EBI 4색 체계 정확히 적용: `>90 → #0053D6`, `70-90 → #65CBF3`, `50-70 → #FFDB13`, `<50 → #FF7D45`
- pLDDT=70 임계선 점선으로 표시
- 500 잔기 초과 시 자동 다운샘플링 + 안내 문구

### STRING 이미지
- URL 패턴: `https://string-db.org/api/image/network?identifiers={gene}&species={taxon}&limit=10&network_flavor=evidence`
- 종(taxon) ID를 organism 이름에서 추론하는 `_guessNCBITaxon()` helper 함수 구현
- 이미지 로드 실패 시 `onerror`로 silently hide

### ChEMBL 구조 이미지
- ADMET 카드에 통합: 기존 테이블 행 → 구조 이미지 포함 카드 레이아웃으로 업그레이드
- ChEMBL ID (`drug_id` 필드)가 "CHEMBL"로 시작할 때만 이미지 표시
- hERG 블로커, CYP3A4 억제제 배지 추가

---

## 향후 과제 (P2, 1주 내)

1. **Human Protein Atlas 조직 발현 히트맵**: `fetchers/protein_atlas.py`에 `fetch_tissue_expression(gene)` 추가 → D3 히트맵 (조직 × 발현 수준)
2. **RCSB 결정 구조 뷰어**: UniProt `dbReferences`에서 PDB ID 추출 → NGL.js로 X-ray 구조 로드. AlphaFold 탭 / 실험 구조 탭 분리

---

## 부록: 기존 V4 회의록 연계

V4에서 결정한 Target Discovery Mode와 V5 UI 업그레이드는 독립적으로 구현된다. Target Discovery 결과 카드(`showTargetDiscoveryReport`)도 V5 디자인 원칙(점수 막대, 배지 체계, 접이식 패널)을 따른다.

---

*회의록 작성: Claude Sonnet 4.6 · 2026-03-20*
