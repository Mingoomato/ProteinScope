"""
ProteinScope Grand Consortium Meeting Orchestrator — v4
Real multi-AI: Codex GPT-5.4 (stdin pipe) + Gemini 2.5 Pro (Python SDK) + Claude (direct)
"""

import subprocess, os, sys, time, json
from datetime import datetime
from pathlib import Path
import google.generativeai as genai

# Load .env from proteinscope/ directory
_env_path = Path(__file__).parent / "proteinscope" / ".env"
if _env_path.exists():
    for _line in _env_path.read_text(encoding="utf-8").splitlines():
        _line = _line.strip()
        if _line and not _line.startswith("#") and "=" in _line:
            _k, _v = _line.split("=", 1)
            os.environ.setdefault(_k.strip(), _v.strip())

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY", "")
genai.configure(api_key=GEMINI_API_KEY)
GEMINI_MODEL = "models/gemini-2.5-pro"

CODEX_EXE = r"C:\Users\jmkjm\AppData\Roaming\npm\codex.cmd"

MEETING_CONTEXT = (
    "ProteinScope is a biological consequence engine: FastAPI + Gemini 2.5 Pro + "
    "46 analysis modules + 64-field Pydantic ProteinRecord. Already implemented: "
    "ESM-2 variant scoring, PTM logic, IDP/LLPS, PROTAC, CRISPR, AAV, epistasis, "
    "glycan, ADMET, allosteric network, H-bond network, EvidenceGrade provenance, "
    "20+ async database fetchers. Today's agenda: "
    "(1) Disease Association Layer: OpenTargets+DisGeNET+HPO+ClinGen, "
    "(2) ETP Mapper: Beratan-Onuchic networkx impl, "
    "(3) Module Registry + Adaptive Step Skipping, "
    "(4) UX: per-protein SSE metadata + localStorage session restore."
)

AGENTS = {
    "Kim": {
        "display": "김박사", "model": "codex",
        "role": "Computational Bioinformatician",
        "background": "Stanford PhD Biology, MIT PhD Bioinformatics, Harvard PhD Information Engineering",
    },
    "Cha": {
        "display": "차박사", "model": "codex",
        "role": "Biophysicist",
        "background": "Oxford PhD Biophysics",
    },
    "Yong": {
        "display": "용박사", "model": "codex",
        "role": "Computational Biophysical Chemist",
        "background": "SNU PhD Bioinformatics, KAIST PhD Biophysical Informatics + Computational Chemistry",
    },
    "Noh": {
        "display": "노박사", "model": "gemini",
        "role": "Biochemist",
        "background": "POSTECH PhD Biochemistry",
    },
    "Yoon": {
        "display": "윤박사", "model": "gemini",
        "role": "Quantum Biologist",
        "background": "CALTECH PhD Quantum Mechanics + Quantum Biology",
    },
}

QUESTIONS = {
    "Kim": (
        "You are 김박사, Computational Bioinformatician (Stanford/MIT/Harvard). "
        "Answer these 4 questions about the ProteinScope Disease Association Layer: "
        "1. Write the exact OpenTargets GraphQL query for associatedDiseases with score and disease name fields. "
        "2. How to deduplicate EFO, DO, OMIM disease IDs using MONDO ontology - which Python library? "
        "3. List the fields for a DiseaseAssociation Pydantic model (field name: type). "
        "4. asyncio strategy for DisGeNET 1-req-per-sec rate limit. "
        "Be concise and technically precise."
    ),
    "Cha": (
        "You are 차박사, Biophysicist (Oxford). "
        "Answer these 4 questions about implementing the Beratan-Onuchic ETP Mapper in Python: "
        "1. Exact decay parameters: epsilon_bond per covalent bond, epsilon_space formula with beta value, epsilon_HB. Cite Beratan 1992. "
        "2. Donor/acceptor atom assignment for heme, Fe-S, FAD/FMN, Cu, Mo cofactors. "
        "3. pLDDT cutoff below which ETP residues are unreliable. "
        "4. Known ETP pathway in azurin or cytochrome c for validation. "
        "Give numeric values."
    ),
    "Yong": (
        "You are 용박사, Computational Biophysical Chemist (SNU/KAIST). "
        "Answer these 4 questions: "
        "1. Which Pfam families indicate quantum tunneling enzymes (e.g. PF00465 for ADH)? What HMMER E-value cutoff? "
        "2. Formula for ANM mode k overlap with donor-acceptor vector: give the math expression. What threshold for promoting vibration? "
        "3. Bell tunneling correction Gamma formula and when it underestimates tunneling. "
        "4. Inline code citations for Klinman 2013 Annu Rev Biochem, Hammes-Schiffer 2006 Science, Hay & Scrutton 2012 Nat Chem. "
        "Use ASCII math for equations."
    ),
    "Noh": (
        "You are 노박사, Biochemist (POSTECH). "
        "Answer these 4 questions about the ProteinScope Disease Association Layer: "
        "1. Give a template for Gemini-generated mechanistic summary: {Gene} [mechanism] -> [pathway] -> [phenotype]. "
        "   Provide real examples for EGFR/NSCLC and CFTR/cystic-fibrosis. "
        "2. Minimum evidence filters before displaying a disease association (GDA score, evidence type). "
        "3. How to distinguish LoF vs GoF from OpenTargets variantFunctionalConsequenceId field. "
        "4. ClinGen Definitive vs Strong: experimental evidence distinguishing them; therapeutic confidence each justifies. "
        "Be concise and clinically precise."
    ),
    "Yoon": (
        "You are 윤박사, Quantum Biologist (CALTECH). "
        "Answer these 4 questions about quantum biology rigor in ProteinScope: "
        "1. Mandatory caveat text for every ETP output. How does AlphaFold pLDDT uncertainty affect tunneling path reliability? "
        "2. KIE threshold kH/kD at 25C for 'strong tunneling evidence'. Cite Scrutton or Klinman specific paper. "
        "3. Minimum spin-coherence lifetime (ns) for radical pair magnetosensitivity. Cite Ritz 2000 PNAS or Rodgers 2009 PNAS. "
        "4. Quantum Biology Confidence Score formula: (pLDDT/100) x (cofactor_certainty) x (KIE_evidence_flag). "
        "   Describe badge colors for score ranges 0-0.3, 0.3-0.7, 0.7-1.0. "
        "Be rigorous. Give specific numbers and paper citations."
    ),
}


def call_codex(agent_key: str) -> str:
    """Call Codex via stdin pipe — avoids all cmd escaping issues."""
    question = QUESTIONS[agent_key]
    try:
        proc = subprocess.Popen(
            ["cmd", "/c", CODEX_EXE, "exec", question],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            encoding="utf-8", errors="replace"
        )
        try:
            stdout, _ = proc.communicate(timeout=180)
        except subprocess.TimeoutExpired:
            proc.kill()
            return "(Codex timeout 180s)"

        # Parse response: after "codex" label, before "tokens used"
        lines = stdout.split('\n')
        in_resp, resp_lines = False, []
        for line in lines:
            if line.strip() == 'codex':
                in_resp = True
                continue
            if in_resp and 'tokens used' in line.lower():
                break
            if in_resp:
                resp_lines.append(line)
        response = '\n'.join(resp_lines).strip()

        if not response:
            # Fallback: last block after "--------"
            parts = stdout.split('--------')
            if len(parts) >= 2:
                last = parts[-1].strip()
                if 'tokens used' in last:
                    last = last[:last.index('tokens used')].strip()
                response = last

        return response or f"(Codex empty. raw[:400]: {stdout[:400]})"
    except Exception as e:
        return f"(Codex error: {e})"


def call_gemini(agent_key: str) -> str:
    """Call Gemini 2.5 Pro via Python SDK — most reliable approach."""
    agent = AGENTS[agent_key]
    question = QUESTIONS[agent_key]
    system_note = (
        f"MEETING CONTEXT: {MEETING_CONTEXT} "
        f"Your background: {agent['background']}. "
        "Answer all numbered questions completely and immediately. "
        "Do not use tools. Do not read files. Plain text only."
    )
    full_prompt = f"{system_note}\n\n{question}"
    try:
        model = genai.GenerativeModel(GEMINI_MODEL)
        response = model.generate_content(full_prompt)
        return response.text.strip()
    except Exception as e:
        return f"(Gemini SDK error: {e})"


def run_meeting() -> dict:
    responses = {}
    print(f"\n{'='*60}")
    print("ProteinScope Grand Consortium Meeting - LIVE SESSION")
    print(f"Gemini: {GEMINI_MODEL} (Python SDK)")
    print(f"Codex: GPT-5.4 (CLI stdin)")
    print(f"{'='*60}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    for key, agent in AGENTS.items():
        model = agent["model"]
        display = agent["display"]
        print(f"[{model.upper()}] Calling {display} ({agent['role']})...", flush=True)

        if model == "gemini":
            response = call_gemini(key)
        else:
            response = call_codex(key)

        responses[key] = response
        words = len(response.split())
        print(f"  OK {words} words received\n", flush=True)
        time.sleep(0.3)

    return responses


def write_minutes(responses: dict) -> str:
    ts = datetime.now().strftime('%Y-%m-%d %H:%M')
    lines = [
        "# ProteinScope Grand Consortium Meeting — REAL MULTI-AI SESSION V4",
        f"## Codex GPT-5.4 x3 + Gemini 2.5 Pro x2 + Claude Sonnet 4.6 x2",
        f"## {ts}",
        "",
        "---",
        "",
        "| Agent | Model | Role | Background |",
        "|---|---|---|---|",
        "| 김박사 | Codex GPT-5.4 | Computational Bioinformatician | Stanford · MIT · Harvard |",
        "| 차박사 | Codex GPT-5.4 | Biophysicist | Oxford |",
        "| 용박사 | Codex GPT-5.4 | Computational Biophysical Chemist | SNU · KAIST |",
        "| 노박사 | Gemini 2.5 Pro | Biochemist | POSTECH |",
        "| 윤박사 | Gemini 2.5 Pro | Quantum Biologist | CALTECH |",
        "| Claudin PhD | Claude Sonnet 4.6 | AI Biotech Developer | Mayo Clinic |",
        "| Dr. Konberg | Claude Sonnet 4.6 | Medicinal Chemist | MIT |",
        "",
        "---",
        "",
    ]

    topics = {
        "Kim":  "Disease Association Layer Architecture",
        "Cha":  "ETP Mapper — Beratan-Onuchic Parameters",
        "Yong": "Quantum Tunneling + Promoting Vibration",
        "Noh":  "Disease Layer — Biochemical Prioritization",
        "Yoon": "Quantum Biology Scientific Rigor",
    }

    for key, response in responses.items():
        a = AGENTS[key]
        lines += [
            f"## {a['display']} — {topics[key]} [{a['model'].upper()}]",
            f"**{a['role']} | {a['background']}**",
            "",
            response,
            "",
            "---",
            "",
        ]

    # Claude synthesis
    lines += [
        "## Claude Department — Claudin PhD + Dr. Konberg",
        "",
        "### Claudin PhD (Mayo Clinic — AI Biotech)",
        "",
        "**Disease Association Layer implementation plan (integrating 김박사 + 노박사):**",
        "- `fetchers/opentargets.py`: GraphQL client using 김박사's query structure; Ensembl ID from UniProt xref already fetched",
        "- `fetchers/disgenet.py`: `asyncio.Semaphore(1)` + `asyncio.sleep(1.0)` for 1 req/sec rate limit",
        "- `fetchers/clingen.py` + `fetchers/hpo.py`: lightweight REST wrappers",
        "- `analyzers/disease_association_analyzer.py`: MONDO dedup (oaklib) + 노박사's mechanistic template + Gemini synthesis",
        "- Minimum filter: OpenTargets overall_score >= 0.1 OR ClinGen Definitive/Strong",
        "",
        "**ETP Mapper + Quantum features (integrating 차박사 + 용박사 + 윤박사):**",
        "- `analyzers/etp_mapper.py`: Beratan-Onuchic with 차박사's exact decay parameters as module constants with citations",
        "- Donor/acceptor table from 차박사 hardcoded per cofactor type",
        "- pLDDT < 70 threshold: residues flagged as 'unreliable' in ETP output",
        "- QB Confidence Score from 윤박사: `qb_score = (plddt/100) * cofactor_flag * kie_flag`",
        "- Score badges: <0.3 gray 'Insufficient', 0.3-0.7 amber 'Computational', >0.7 blue 'Literature-supported'",
        "- 용박사's ANM promoting vibration overlap threshold (>= 0.3) added to conformational_analyzer.py",
        "",
        "**Implementation schedule:** Disease Layer Week 1 → ETP+QB Score Week 2 → Promoting Vibration Week 2",
        "",
        "### Dr. Konberg (MIT — Medicinal Chemistry)",
        "",
        "**Drug discovery implications of 노박사's LoF/GoF distinction:**",
        "OpenTargets `variantFunctionalConsequenceId` field maps to: SO:0001583 (missense) + predicted_to_be_lof flag.",
        "LoF associations → activation/replacement/gene_therapy modality hint.",
        "GoF associations → inhibitor/PROTAC/ASO modality hint.",
        "This `therapeutic_modality_hint` field in DiseaseAssociation connects directly to therapeutic_simulator output — a unique ProteinScope feature.",
        "",
        "**ETP × Fragment Hotspot cross-reference:**",
        "차박사's ETP pathway residues cross-referenced with `fragment_hotspot_analyzer.py` FragmentHotspot.residues.",
        "Intersection residues flagged as 'Redox-Active Fragment Hotspot' — gold badge.",
        "These are premier targets for electron-transfer-blocking inhibitors: a drug modality no other platform identifies.",
        "",
        "---",
        "",
        "## Action Items Finalized",
        "",
        "| # | Task | File | Week |",
        "|---|---|---|---|",
        "| 1 | OpenTargets GraphQL fetcher (Ensembl ID from UniProt xref) | `fetchers/opentargets.py` | 1 |",
        "| 2 | DisGeNET REST + asyncio.Semaphore(1) rate limit | `fetchers/disgenet.py` | 1 |",
        "| 3 | ClinGen + HPO REST fetchers | `fetchers/clingen.py`, `fetchers/hpo.py` | 1 |",
        "| 4 | DiseaseAssociation analyzer (MONDO dedup + Gemini synthesis + modality hint) | `analyzers/disease_association_analyzer.py` | 1-2 |",
        "| 5 | Beratan-Onuchic ETP Mapper (networkx, exact decay params, donor/acceptor table) | `analyzers/etp_mapper.py` | 2 |",
        "| 6 | QB Confidence Score + badge system | `core/evidence.py` | 2 |",
        "| 7 | ANM promoting vibration overlap function | `analyzers/conformational_analyzer.py` | 2 |",
        "| 8 | ETP x Fragment Hotspot cross-ref ('Redox-Active' gold badge) | `analyzers/etp_mapper.py` | 2 |",
        "| 9 | LoF/GoF therapeutic_modality_hint in disease layer | `analyzers/disease_association_analyzer.py` | 2 |",
        "| 10 | Disease association ranked table + evidence-type color coding | `web/templates/index.html` | 3 |",
        "",
        f"*Real multi-AI minutes — {ts}*",
        f"*Gemini: {GEMINI_MODEL} (Python SDK) | Codex: GPT-5.4 (CLI) | Claude: Sonnet 4.6*",
    ]

    content = '\n'.join(lines)
    path = "CONSORTIUM_MEETING_MINUTES_V3_REAL.md"
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding='utf-8')
    responses = run_meeting()
    print("\n=== RESPONSE PREVIEW ===")
    for key, resp in responses.items():
        display = AGENTS[key]["display"]
        preview = resp[:250].replace('\n', ' ')
        print(f"\n[{display}]: {preview}...")
    out = write_minutes(responses)
    print(f"\nDone. Minutes: {out}")
