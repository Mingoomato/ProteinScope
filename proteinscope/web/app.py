"""ProteinScope Web UI — FastAPI application.

Run with:
    python main.py serve
    or: uvicorn web.app:app --reload --port 8000
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
import uuid
from pathlib import Path
from typing import AsyncGenerator

from fastapi import BackgroundTasks, FastAPI, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from dotenv import load_dotenv

# Ensure proteinscope package root is on path
_ROOT = str(Path(__file__).parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

load_dotenv(Path(_ROOT) / ".env")

app = FastAPI(title="ProteinScope", version="1.0.0")

_TEMPLATES = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))
_STATIC = Path(__file__).parent / "static"

app.mount("/static", StaticFiles(directory=str(_STATIC)), name="static")

# Serve generated reports (pathway images, etc.)
_REPORTS_DIR = Path(_ROOT) / "reports"
_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
app.mount("/reports", StaticFiles(directory=str(_REPORTS_DIR)), name="reports")

# Serve pathway diagram images saved by diagram_fetcher (save_dir="./diagrams")
_DIAGRAMS_DIR = Path(_ROOT) / "diagrams"
_DIAGRAMS_DIR.mkdir(parents=True, exist_ok=True)
app.mount("/diagrams", StaticFiles(directory=str(_DIAGRAMS_DIR)), name="diagrams")

# Mount API v1 router
from api.v1.router import router as _api_v1  # noqa: E402
app.include_router(_api_v1, prefix="/api/v1")

# ---------------------------------------------------------------------------
# In-memory job store
# ---------------------------------------------------------------------------

_jobs: dict[str, dict] = {}          # job_id -> metadata
_queues: dict[str, asyncio.Queue] = {}  # job_id -> SSE queue


def _new_job() -> str:
    jid = str(uuid.uuid4())
    _jobs[jid] = {"status": "pending", "accession": None, "report_html": None,
                  "report_md": None, "report_pdf": None, "error": None}
    _queues[jid] = asyncio.Queue()
    return jid


async def _push(jid: str, msg: dict):
    await _queues[jid].put(msg)


# ---------------------------------------------------------------------------
# Background query runner
# ---------------------------------------------------------------------------

async def _run_query(
    jid: str,
    query: str,
    organism: str | None,
    cross_species: bool,
    include_pathways: bool,
    ai_summary: bool,
    no_images: bool,
    min_pgx_evidence: str,
    interaction_score: float,
    session_id: str = "",
):
    from core.query_engine import _resolve_accession, fetch_all, ProteinSuggestionsError, TaxonAmbiguityError, SpeciesFallbackError, OffTopicQueryError, ReverseGeneticsRequest, RNAiRequest, PathwayAnalysisRequest
    from core.cache import Cache
    from output.markdown_writer import write_markdown
    from output.pdf_writer import write_pdf

    try:
        await _push(jid, {"type": "step", "step": 1, "message": "Resolving protein identifier…"})
        cache = Cache()

        # Handle compatibility redirect — raised when a query implicitly
        # targets a cross-species comparison rather than a single protein.
        try:
            accession = await _resolve_accession(query, organism)
        except Exception as exc:
            # Check whether the exception signals a compatibility request
            exc_type_name = type(exc).__name__
            if exc_type_name == "ProteinCompatibilityRequest":
                _jobs[jid]["status"] = "compatibility_redirect"
                payload = {
                    "type": "compatibility_redirect",
                    "source_gene": getattr(exc, "source_gene", query),
                    "source_organism": getattr(exc, "source_organism", organism),
                    "target_organism": getattr(exc, "target_organism", None),
                    "target_cell": getattr(exc, "target_cell", None),
                }
                _jobs[jid].update(payload)
                await _push(jid, payload)
                return
            if exc_type_name == "ReverseGeneticsRequest":
                _jobs[jid]["status"] = "reverse_genetics_redirect"
                payload = {
                    "type": "reverse_genetics_redirect",
                    "gene": getattr(exc, "gene", query),
                    "mutation_type": getattr(exc, "mutation_type", "general"),
                    "specific_variant": getattr(exc, "specific_variant", None),
                    "organism": getattr(exc, "organism", "Homo sapiens"),
                }
                _jobs[jid].update(payload)
                await _push(jid, payload)
                return
            if exc_type_name == "PathwayAnalysisRequest":
                _jobs[jid]["status"] = "pathway_plan"
                payload = {
                    "type": "pathway_plan",
                    "query": query,
                    "plan": getattr(exc, "plan", {}),
                }
                _jobs[jid].update(payload)
                await _push(jid, payload)
                return
            if exc_type_name == "RNAiRequest":
                _jobs[jid]["status"] = "rnai_redirect"
                payload = {
                    "type": "rnai_redirect",
                    "gene": getattr(exc, "gene", query),
                    "rnai_type": getattr(exc, "rnai_type", "siRNA"),
                    "cell_type": getattr(exc, "cell_type", None),
                    "organism": getattr(exc, "organism", "Homo sapiens"),
                }
                _jobs[jid].update(payload)
                await _push(jid, payload)
                return
            raise

        _jobs[jid]["accession"] = accession
        await _push(jid, {"type": "step", "step": 2,
                          "message": f"Resolved → {accession}. Fetching protein data…"})

        record = await fetch_all(
            accession, organism, cross_species, cache,
            include_pathways=include_pathways,
            no_images=no_images,
            min_pgx_evidence=min_pgx_evidence,
            interaction_score=interaction_score,
            ai_summary=ai_summary,
        )
        cache.close()

        await _push(jid, {"type": "step", "step": 3, "message": "Rendering report…"})

        out_dir = Path(_ROOT) / "reports" / jid
        out_dir.mkdir(parents=True, exist_ok=True)
        out_base = str(out_dir / "report")

        md_path = write_markdown(record, out_base)
        pdf_path = write_pdf(record, out_base)

        # Convert markdown to HTML for inline display
        import markdown as md_lib, re as _re
        md_text = md_path.read_text(encoding="utf-8")

        # markdown_writer.p() adds a trailing \n per line; "\n".join() then
        # inserts another \n between elements, producing blank lines between
        # table rows.  The `tables` extension requires contiguous rows, so
        # remove the spurious blank line between any two pipe-starting lines.
        md_text = _re.sub(r'(\|[^\n]*)\n\n(?=\|)', r'\1\n', md_text)

        report_html = md_lib.markdown(
            md_text,
            extensions=["tables", "fenced_code", "sane_lists"],
        )

        # Replace absolute local image paths (from diagram_image_path) with
        # web-accessible URLs so <img> tags render.
        def _fix_img_src(m: "_re.Match") -> str:
            raw = m.group(1)
            # Normalise Windows backslashes
            norm = raw.replace("\\", "/")
            # Case 1: path under reports/{jid}/
            marker = f"/reports/{jid}/"
            idx = norm.find(marker)
            if idx != -1:
                rel = norm[idx + len(marker):]
                return f'src="/reports/{jid}/{rel}"'
            # Case 2: path under diagrams/ (saved by diagram_fetcher)
            diag_marker = "diagrams/"
            idx2 = norm.find(diag_marker)
            if idx2 != -1:
                rel2 = norm[idx2 + len(diag_marker):]
                return f'src="/diagrams/{rel2}"'
            return m.group(0)  # unchanged

        report_html = _re.sub(r'src="([^"]+)"', _fix_img_src, report_html)

        # Serialise structured data for the frontend (evidence badges, ESM-2 bars, network)
        _jobs[jid].update({
            "status": "done",
            "report_html": report_html,
            "report_md": str(md_path),
            "report_pdf": str(pdf_path),
            "protein_name": record.protein_name,
            "gene_name": record.gene_name,
            "organism": record.organism,
            "uniprot_id": record.uniprot_id,
            # Structured data used by frontend JS (evidence badges, ESM-2, PPI network)
            "protein_interactions": [i.model_dump() for i in record.protein_interactions],
            "clinical_variants": [c.model_dump() for c in record.clinical_variants],
            "drug_interactions": [d.model_dump() for d in record.drug_interactions],
            "variant_fitness_scores": [v.model_dump() for v in record.variant_fitness_scores],
            "alphafold_plddt_score": record.alphafold_plddt_score,
            "alphafold_pdb_url": record.alphafold_pdb_url,
            "canonical_sequence": record.canonical_sequence,  # needed by /api/esm_score
            # P2: PTM Logic
            "ptm_logic": record.ptm_logic.model_dump() if record.ptm_logic else None,
            # P3: IDP/LLPS
            "idp_analysis": record.idp_analysis.model_dump() if record.idp_analysis else None,
            # P4: Metalloenzyme
            "metalloenzyme_analysis": record.metalloenzyme_analysis.model_dump() if record.metalloenzyme_analysis else None,
            # P5: Expression Compatibility
            "expression_compatibility": record.expression_compatibility.model_dump() if record.expression_compatibility else None,
            # P6: Druggability
            "druggability_report": record.druggability_report.model_dump() if record.druggability_report else None,
            # Ph2: Conformational
            "conformational_analysis": record.conformational_analysis.model_dump() if record.conformational_analysis else None,
            # P9: Inverse folding designs
            "inverse_designs": [d.model_dump() for d in record.inverse_designs],
            # P10: Therapeutic Decision Simulator
            "therapeutic_decision": record.therapeutic_decision.model_dump() if record.therapeutic_decision else None,
            # P11: Antigen Discovery (per-protein HPA expression profile)
            "antigen_discovery": record.antigen_discovery.model_dump() if record.antigen_discovery else None,
            # Layer 1: Molecular State
            "allosteric_network": record.allosteric_network.model_dump() if record.allosteric_network else None,
            "epistasis_report": record.epistasis_report.model_dump() if record.epistasis_report else None,
            "binding_energy_estimates": [e.model_dump() for e in record.binding_energy_estimates],
            "observables_prediction": record.observables_prediction.model_dump() if record.observables_prediction else None,
            # Layer 2: Predictive Mechanism
            "active_learning_recommendation": record.active_learning_recommendation.model_dump() if record.active_learning_recommendation else None,
            # Layer 3: Engineering Feasibility
            "glycan_analysis": record.glycan_analysis.model_dump() if record.glycan_analysis else None,
            "admet_report": record.admet_report.model_dump() if record.admet_report else None,
            "directed_evolution_plan": record.directed_evolution_plan.model_dump() if record.directed_evolution_plan else None,
            "protac_feasibility": record.protac_feasibility.model_dump() if record.protac_feasibility else None,
            "selectivity_landscape": record.selectivity_landscape.model_dump() if record.selectivity_landscape else None,
            "fragment_hotspot_map": record.fragment_hotspot_map.model_dump() if record.fragment_hotspot_map else None,
            "aav_design": record.aav_design.model_dump() if record.aav_design else None,
            "crispr_plan": record.crispr_plan.model_dump() if record.crispr_plan else None,
            "genetic_circuit_design": record.genetic_circuit_design.model_dump() if record.genetic_circuit_design else None,
            # Layer 3+: Additional Engineering
            "covalent_inhibitor_design": record.covalent_inhibitor_design.model_dump() if record.covalent_inhibitor_design else None,
            "engineering_strategy": record.engineering_strategy.model_dump() if record.engineering_strategy else None,
            "hbond_network": record.hbond_network.model_dump() if record.hbond_network else None,
        })
        await _push(jid, {"type": "done", "message": "Report ready!"})

        # P8: index protein for concept search
        try:
            from analyzers.embedding_search import index_protein
            index_protein(
                gene=record.gene_name,
                accession=record.uniprot_id,
                organism=record.organism,
                description=record.function_description[:500],
                go_terms=record.biological_process[:10] + record.molecular_function[:10],
            )
        except Exception:
            pass

        # Ph4: update session history
        if session_id:
            try:
                from core.session_manager import session_manager as _sm
                _jobs[jid]["job_id"] = jid
                _sm.update(session_id, _jobs[jid])
            except Exception:
                pass

    except ProteinSuggestionsError as exc:
        _jobs[jid]["status"] = "suggestions"
        await _push(jid, {
            "type": "suggestions",
            "query": exc.query,
            "proteins": exc.suggestions,
            "organism": exc.organism,
        })

    except TaxonAmbiguityError as e:
        await _push(jid, {
            "type": "taxon_ambiguity",
            "organism": e.organism,
            "taxon_rank": e.taxon_rank,
            "species": e.species,
        })

    except SpeciesFallbackError as e:
        _jobs[jid]["status"] = "species_fallback"
        _jobs[jid]["species_fallback_accession"] = e.accession
        await _push(jid, {
            "type": "species_fallback",
            "gene": e.gene,
            "requested_organism": e.requested_organism,
            "actual_organism": e.actual_organism,
            "accession": e.accession,
        })

    except OffTopicQueryError as e:
        _jobs[jid]["status"] = "off_topic"
        await _push(jid, {
            "type": "off_topic",
            "query": e.query,
            "reason": e.reason,
        })

    except Exception as exc:
        _jobs[jid]["status"] = "error"
        _jobs[jid]["error"] = str(exc)
        await _push(jid, {"type": "error", "message": str(exc)})


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.get("/health")
async def health():
    """Health check endpoint — returns status plus GPU device info."""
    try:
        from utils.device import device_info
        gpu = device_info()
    except Exception:
        gpu = {"device": "unknown"}
    return {"status": "ok", "version": "1.0.0", **gpu}


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _TEMPLATES.TemplateResponse("index.html", {"request": request})


@app.post("/api/search")
async def start_search(request: Request, background_tasks: BackgroundTasks):
    form = await request.form()
    query = (form.get("query") or "").strip()
    if not query:
        return {"error": "Query is required"}

    organism = form.get("organism") or None
    cross_species = form.get("cross_species") == "true"
    include_pathways = form.get("no_pathways") != "true"
    ai_summary = form.get("ai_summary") == "true"
    no_images = form.get("no_images") == "true"
    min_pgx_evidence = form.get("min_pgx_evidence") or "4"
    interaction_score = float(form.get("interaction_score") or "0.4")

    # Ph4: Session context — read or generate session ID
    session_id = request.headers.get("X-Session-ID", "").strip() or str(uuid.uuid4())

    jid = _new_job()
    background_tasks.add_task(
        _run_query, jid, query, organism, cross_species,
        include_pathways, ai_summary, no_images,
        min_pgx_evidence, interaction_score, session_id,
    )
    return JSONResponse({"job_id": jid, "session_id": session_id},
                        headers={"X-Session-ID": session_id})


@app.get("/api/progress/{job_id}")
async def progress_stream(job_id: str):
    async def event_generator() -> AsyncGenerator[str, None]:
        queue = _queues.get(job_id)
        if not queue:
            yield f"data: {json.dumps({'type': 'error', 'message': 'Job not found'})}\n\n"
            return
        while True:
            try:
                msg = await asyncio.wait_for(queue.get(), timeout=120)
                yield f"data: {json.dumps(msg)}\n\n"
                if msg.get("type") in ("done", "error", "cancelled", "compatibility_redirect", "compatibility_result", "taxon_ambiguity", "suggestions", "species_fallback", "off_topic", "reverse_genetics_redirect", "rnai_redirect", "reverse_genetics_result", "rnai_result", "pathway_plan", "pathway_summary_result", "esm_result"):
                    break
            except asyncio.TimeoutError:
                yield "data: {\"type\": \"ping\"}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no"})


@app.post("/api/cancel/{job_id}")
async def cancel_job(job_id: str):
    """Mark a job as cancelled so its background task stops sending events."""
    job = _jobs.get(job_id)
    if job:
        job["status"] = "cancelled"
        # Push a sentinel so the SSE generator loop exits cleanly if still open
        queue = _queues.get(job_id)
        if queue:
            await queue.put({"type": "cancelled"})
    return {"ok": True}


@app.get("/api/report/{job_id}")
async def get_report(job_id: str):
    job = _jobs.get(job_id)
    if not job:
        return {"error": "Job not found"}
    return job


@app.post("/api/esm_score")
async def start_esm_score(request: Request, background_tasks: BackgroundTasks):
    """Launch an ESM-2 fitness landscape scan for a protein.

    Form fields:
        job_id     — the completed protein job whose canonical sequence to use
        positions  — optional comma-separated 1-based positions to scan
                     (default: all positions — slow for long sequences)
        model_name — optional ESM-2 model name (default: esm2_t12_35M_UR50D for speed)

    Returns: {esm_job_id}
    The SSE stream emits {"type": "esm_result", "scores": [...VariantScore dicts...]}
    when complete.
    """
    form = await request.form()
    job_id = (form.get("job_id") or "").strip()
    positions_str = (form.get("positions") or "").strip()
    model_name = (form.get("model_name") or "esm2_t12_35M_UR50D").strip()

    source_job = _jobs.get(job_id)
    if not source_job:
        return JSONResponse({"error": "Source job not found"}, status_code=404)

    # Extract canonical sequence from the job's clinical data or from ProteinRecord
    # The sequence isn't stored in _jobs yet — we re-fetch it from the report_md metadata.
    # Simpler: return error so the frontend knows to call /api/esm_score only with valid job.
    # We store canonical_sequence in a later step; for now use a lazy approach.

    esm_jid = _new_job()
    background_tasks.add_task(
        _run_esm_score, esm_jid, source_job, positions_str, model_name
    )
    return {"esm_job_id": esm_jid}


async def _run_esm_score(esm_jid: str, source_job: dict, positions_str: str, model_name: str):
    """Background task: run score_all_possible_variants and push SSE result."""
    try:
        canonical_sequence = source_job.get("canonical_sequence", "")
        if not canonical_sequence:
            await _push(esm_jid, {"type": "error", "message": "No sequence available — run a protein search first"})
            return

        await _push(esm_jid, {"type": "step", "message": f"Loading ESM-2 model ({model_name})…"})

        from analyzers.variant_scorer import score_all_possible_variants
        import asyncio

        positions = None
        if positions_str:
            try:
                positions = [int(p.strip()) for p in positions_str.split(",") if p.strip()]
            except ValueError:
                pass

        loop = asyncio.get_event_loop()
        scores = await loop.run_in_executor(
            None, score_all_possible_variants, canonical_sequence, positions, model_name
        )

        score_dicts = [s.model_dump() for s in scores]
        _jobs[esm_jid]["esm_scores"] = score_dicts
        await _push(esm_jid, {"type": "esm_result", "scores": score_dicts})

    except Exception as exc:
        await _push(esm_jid, {"type": "error", "message": f"ESM-2 scoring failed: {exc}"})


@app.get("/api/download/{job_id}/pdf")
async def download_pdf(job_id: str):
    job = _jobs.get(job_id)
    if not job or not job.get("report_pdf"):
        return HTMLResponse("Not found", status_code=404)
    fname = f"{job.get('gene_name', 'report')}_ProteinScope.pdf"
    return FileResponse(job["report_pdf"], media_type="application/pdf",
                        filename=fname)


@app.get("/api/download/{job_id}/markdown")
async def download_markdown(job_id: str):
    job = _jobs.get(job_id)
    if not job or not job.get("report_md"):
        return HTMLResponse("Not found", status_code=404)
    fname = f"{job.get('gene_name', 'report')}_ProteinScope.md"
    return FileResponse(job["report_md"], media_type="text/markdown",
                        filename=fname)


@app.post("/api/compatibility")
async def start_compatibility(request: Request, background_tasks: BackgroundTasks):
    """Submit a cross-species compatibility analysis job.

    Form fields: source_gene, source_organism, target_organism, target_cell (optional).
    Returns: {job_id}
    """
    form = await request.form()
    source_gene = (form.get("source_gene") or "").strip()
    source_organism = (form.get("source_organism") or "").strip()
    target_organism = (form.get("target_organism") or "").strip()
    target_cell = (form.get("target_cell") or None)

    if not source_gene or not source_organism or not target_organism:
        return JSONResponse(
            {"error": "source_gene, source_organism, and target_organism are required"},
            status_code=422,
        )

    jid = _new_job()
    _jobs[jid]["source_gene"] = source_gene
    _jobs[jid]["source_organism"] = source_organism
    _jobs[jid]["target_organism"] = target_organism
    _jobs[jid]["target_cell"] = target_cell

    async def _push_step(msg: str):
        await _push(jid, {"type": "step", "message": msg})

    async def _run_compat(job_id: str):
        try:
            from core.query_engine import run_compatibility_query  # noqa: WPS433
            result = await run_compatibility_query(
                source_gene=source_gene,
                source_organism=source_organism,
                target_organism=target_organism,
                target_cell=target_cell,
                step_cb=_push_step,
            )
            # result is now a dict {"compatibility": ..., "insertion_impact": ...}
            result_data = result if isinstance(result, dict) else (
                result.model_dump() if hasattr(result, "model_dump") else vars(result)
            )
            _jobs[job_id].update({
                "status": "done",
                "compatibility_result": result_data,
            })
            await _push(job_id, {"type": "compatibility_result", "report": result_data})
        except ImportError:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = (
                "run_compatibility_query is not yet implemented in core.query_engine"
            )
            await _push(job_id, {
                "type": "error",
                "message": _jobs[job_id]["error"],
            })
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)
            await _push(job_id, {"type": "error", "message": str(exc)})

    background_tasks.add_task(_run_compat, jid)
    return {"job_id": jid}


@app.post("/api/reverse_genetics")
async def start_reverse_genetics(request: Request, background_tasks: BackgroundTasks):
    """Submit a reverse genetics analysis job.

    Form fields: gene, mutation_type, organism (optional), specific_variant (optional).
    Returns: {job_id}
    """
    form = await request.form()
    gene = (form.get("gene") or "").strip()
    mutation_type = (form.get("mutation_type") or "general").strip()
    organism = (form.get("organism") or "Homo sapiens").strip()
    specific_variant = form.get("specific_variant") or None

    if not gene:
        return JSONResponse({"error": "gene is required"}, status_code=422)

    jid = _new_job()
    _jobs[jid]["gene"] = gene
    _jobs[jid]["mutation_type"] = mutation_type
    _jobs[jid]["organism"] = organism

    async def _push_step(msg: str):
        await _push(jid, {"type": "step", "message": msg})

    async def _run_rg(job_id: str):
        try:
            from core.query_engine import run_reverse_genetics_query
            result = await run_reverse_genetics_query(
                gene=gene,
                mutation_type=mutation_type,
                organism=organism,
                specific_variant=specific_variant,
                step_cb=_push_step,
            )
            _jobs[job_id].update({"status": "done", "reverse_genetics_result": result})
            await _push(job_id, {"type": "reverse_genetics_result", "report": result})
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)
            await _push(job_id, {"type": "error", "message": str(exc)})

    background_tasks.add_task(_run_rg, jid)
    return {"job_id": jid}


@app.post("/api/rnai")
async def start_rnai(request: Request, background_tasks: BackgroundTasks):
    """Submit an RNAi experiment design job.

    Form fields: gene, rnai_type (siRNA|shRNA|miRNA|general), organism (optional),
                 cell_type (optional).
    Returns: {job_id}
    """
    form = await request.form()
    gene = (form.get("gene") or "").strip()
    rnai_type = (form.get("rnai_type") or "siRNA").strip()
    organism = (form.get("organism") or "Homo sapiens").strip()
    cell_type = form.get("cell_type") or None

    if not gene:
        return JSONResponse({"error": "gene is required"}, status_code=422)

    jid = _new_job()
    _jobs[jid]["gene"] = gene
    _jobs[jid]["rnai_type"] = rnai_type
    _jobs[jid]["organism"] = organism

    async def _push_step(msg: str):
        await _push(jid, {"type": "step", "message": msg})

    async def _run_rnai(job_id: str):
        try:
            from core.query_engine import run_rnai_query
            result = await run_rnai_query(
                gene=gene,
                rnai_type=rnai_type,
                organism=organism,
                cell_type=cell_type,
                step_cb=_push_step,
            )
            _jobs[job_id].update({"status": "done", "rnai_result": result})
            await _push(job_id, {"type": "rnai_result", "report": result})
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)
            await _push(job_id, {"type": "error", "message": str(exc)})

    background_tasks.add_task(_run_rnai, jid)
    return {"job_id": jid}


@app.post("/api/pathway_summary")
async def pathway_summary(request: Request, background_tasks: BackgroundTasks):
    """Generate a systems-biology synthesis of all completed pathway engineering steps.

    JSON body: {plan: {...}, completed_steps: [{id, title, result}, ...]}
    Returns: {job_id}
    """
    body = await request.json()
    plan = body.get("plan", {})
    completed_steps = body.get("completed_steps", [])

    jid = _new_job()

    async def _run_summary(job_id: str):
        try:
            from core.gemini_interpreter import generate_pathway_summary
            summary_text = await generate_pathway_summary(plan, completed_steps)
            _jobs[job_id].update({"status": "done", "pathway_summary": summary_text})
            await _push(job_id, {"type": "pathway_summary_result", "summary": summary_text})
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)
            await _push(job_id, {"type": "error", "message": str(exc)})

    background_tasks.add_task(_run_summary, jid)
    return {"job_id": jid}


@app.post("/api/export/pdf")
async def export_pdf(request: Request):
    """Convert a markdown report to PDF and return as a download.

    JSON body: {title: str, markdown: str}
    """
    from fastapi.responses import Response as _Response
    body = await request.json()
    title = body.get("title", "ProteinScope Report")
    markdown_text = body.get("markdown", "")
    from output.pdf_writer import generate_markdown_pdf
    pdf_bytes = generate_markdown_pdf(title, markdown_text)
    safe_name = title.replace(" ", "_").replace("/", "-")[:80]
    return _Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{safe_name}.pdf"'},
    )


@app.post("/api/concept_search")
async def concept_search(request: Request):
    """Search the protein index by natural language concept.

    Form fields:
        query  — the concept query (e.g. "disordered kinase with nuclear localization")
        top_k  — number of results (default 10)

    Returns: {results: [{gene, accession, organism, description, score}], index_size: int}
    """
    form = await request.form()
    query = (form.get("query") or "").strip()
    top_k = int(form.get("top_k") or "10")
    if not query:
        return JSONResponse({"error": "Query is required"}, status_code=400)
    try:
        from analyzers.embedding_search import search_by_concept, get_index_size
        results = search_by_concept(query, top_k=top_k)
        return {"results": results, "index_size": get_index_size()}
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)


@app.post("/api/inverse_fold")
async def start_inverse_fold(request: Request, background_tasks: BackgroundTasks):
    """Submit an inverse-folding design job (ESM-IF1).

    Form fields: job_id (required — must be a completed protein query job),
                 n_designs (default 5), target_host (default "E.coli")
    Returns: {job_id}  — poll /api/progress/{job_id} for inverse_fold_result
    """
    form = await request.form()
    source_job_id = (form.get("job_id") or "").strip()
    n_designs = int(form.get("n_designs") or "5")
    target_host = (form.get("target_host") or "E.coli").strip()

    if not source_job_id or source_job_id not in _jobs:
        return JSONResponse({"error": "Valid job_id required"}, status_code=422)

    source = _jobs[source_job_id]
    af_pdb_url = source.get("alphafold_pdb_url") or ""
    af_plddt = source.get("alphafold_plddt_score")
    gene = source.get("gene_name", "unknown")

    jid = _new_job()
    _jobs[jid]["gene"] = gene
    _jobs[jid]["source_job_id"] = source_job_id

    async def _step(msg: str):
        await _push(jid, {"type": "step", "message": msg})

    async def _run_if(job_id: str):
        try:
            from analyzers.inverse_folder import design_sequences
            designs = await design_sequences(
                gene=gene,
                af_pdb_url=af_pdb_url,
                n_designs=n_designs,
                target_host=target_host,
                af_plddt=af_plddt,
                step_cb=_step,
            )
            result = [d.model_dump() for d in designs]
            _jobs[job_id].update({"status": "done", "inverse_fold_result": result})
            await _push(job_id, {"type": "inverse_fold_result", "designs": result, "gene": gene})
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)
            await _push(job_id, {"type": "error", "message": str(exc)})

    background_tasks.add_task(_run_if, jid)
    return {"job_id": jid}


@app.post("/api/alignment")
async def compute_alignment(request: Request):
    """Compute pairwise sequence alignment and return gap-aligned strings.

    JSON body: {seq_a: str, seq_b: str, seq_type?: "auto"|"protein"|"nucleotide"}
    seq_type defaults to "auto" — auto-detects DNA/RNA vs protein.
    Response: {aligned_query: str, aligned_target: str, identity_pct: float,
               similarity_pct: float, seq_type: str}
    """
    body = await request.json()
    seq_a = (body.get("seq_a") or "").strip()
    seq_b = (body.get("seq_b") or "").strip()
    seq_type = (body.get("seq_type") or "auto").strip().lower()

    if not seq_a or not seq_b:
        return JSONResponse({"error": "seq_a and seq_b are required"}, status_code=422)
    if len(seq_a) < 5 or len(seq_b) < 5:
        return JSONResponse({"error": "Sequences must be at least 5 residues long"}, status_code=422)

    try:
        from analyzers.sequence_aligner import (
            align_pairwise,
            align_pairwise_nucleotide,
            align_pairwise_auto,
            _is_nucleotide_sequence,
        )
        if seq_type == "nucleotide":
            result = align_pairwise_nucleotide(seq_a, seq_b)
            detected_type = "nucleotide"
        elif seq_type == "protein":
            result = align_pairwise(seq_a, seq_b)
            detected_type = "protein"
        else:  # auto
            is_nucl = _is_nucleotide_sequence(seq_a) and _is_nucleotide_sequence(seq_b)
            result = align_pairwise_nucleotide(seq_a, seq_b) if is_nucl else align_pairwise(seq_a, seq_b)
            detected_type = "nucleotide" if is_nucl else "protein"

        if result is None:
            return JSONResponse({"error": "Alignment failed — check sequence format"}, status_code=422)
        return {
            "aligned_query": result.aligned_query,
            "aligned_target": result.aligned_target,
            "identity_pct": result.identity_pct,
            "similarity_pct": result.similarity_pct,
            "seq_type": detected_type,
        }
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)


@app.post("/api/antigen_discovery")
async def start_antigen_discovery(request: Request):
    """Screen Human Protein Atlas for tumor antigen candidates (P11).

    Form fields:
        tumor_type (str, required)  — e.g. "melanoma", "glioblastoma"
        modality   (str, optional)  — CAR_T | antibody | ADC | TCE (default CAR_T)

    Returns:
        AntigenDiscoveryReport JSON with scored candidates and Gemini rationale.
    """
    form = await request.form()
    tumor_type = (form.get("tumor_type") or "").strip()
    modality = (form.get("modality") or "CAR_T").strip()

    if not tumor_type:
        return JSONResponse({"error": "tumor_type is required"}, status_code=422)

    try:
        from fetchers.protein_atlas import fetch_hpa_tumor_expression
        from analyzers.antigen_discovery_analyzer import run_antigen_discovery
        hpa_data = await fetch_hpa_tumor_expression(tumor_type)
        report = await run_antigen_discovery(
            tumor_type=tumor_type,
            modality=modality,
            hpa_data=hpa_data,
        )
        return JSONResponse(report.model_dump())
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)
