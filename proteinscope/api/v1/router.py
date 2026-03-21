"""REST API v1 — programmatic access to ProteinScope.

Mount this router in web/app.py:
    from api.v1.router import router as api_v1
    app.include_router(api_v1, prefix="/api/v1")

Endpoints:
    POST /api/v1/analyze          Submit single protein analysis
    POST /api/v1/compare          Submit compatibility analysis
    GET  /api/v1/jobs/{job_id}    Poll job status
    GET  /api/v1/jobs/{job_id}/report  Get JSON report
    GET  /health                  Health check (also mounted at top level)
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel

router = APIRouter()

# ---------------------------------------------------------------------------
# Request / response schemas
# ---------------------------------------------------------------------------


class AnalyzeRequest(BaseModel):
    query: str
    organism: Optional[str] = None
    ai_summary: bool = True
    no_pathways: bool = False


class CompareRequest(BaseModel):
    source_gene: str
    source_organism: str
    target_organism: str
    target_cell: Optional[str] = None


class JobCreatedResponse(BaseModel):
    job_id: str


class JobStatusResponse(BaseModel):
    job_id: str
    status: str
    created_at: Optional[str] = None
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Helpers — share the _jobs store from web.app
# ---------------------------------------------------------------------------


def _get_jobs() -> dict:
    """Import _jobs lazily to avoid circular-import at module load time."""
    from web.app import _jobs, _new_job  # noqa: WPS433
    return _jobs, _new_job


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@router.get("/health", tags=["meta"])
async def health():
    """Health check — also reachable at GET /health via the main app."""
    return {"status": "ok", "version": "1.0.0"}


@router.post("/analyze", response_model=JobCreatedResponse, tags=["jobs"])
async def analyze(body: AnalyzeRequest, background_tasks: BackgroundTasks):
    """Submit a single-protein analysis job.

    Returns a ``job_id`` that can be polled via ``GET /api/v1/jobs/{job_id}``.
    """
    from web.app import _jobs, _new_job, _run_query  # noqa: WPS433

    if not body.query.strip():
        raise HTTPException(status_code=422, detail="query must not be blank")

    jid = _new_job()
    _jobs[jid]["created_at"] = datetime.now(timezone.utc).isoformat()

    background_tasks.add_task(
        _run_query,
        jid,
        body.query.strip(),
        body.organism,
        False,                    # cross_species
        not body.no_pathways,     # include_pathways
        body.ai_summary,
        False,                    # no_images
        "4",                      # min_pgx_evidence
        0.4,                      # interaction_score
    )
    return JobCreatedResponse(job_id=jid)


@router.post("/compare", response_model=JobCreatedResponse, tags=["jobs"])
async def compare(body: CompareRequest, background_tasks: BackgroundTasks):
    """Submit a cross-species compatibility analysis job.

    ``run_compatibility_query`` will be called when it exists in
    ``core.query_engine``; until then the job immediately errors with a
    descriptive message so callers receive useful feedback.
    """
    from web.app import _jobs, _new_job  # noqa: WPS433

    jid = _new_job()
    _jobs[jid]["created_at"] = datetime.now(timezone.utc).isoformat()
    _jobs[jid]["source_gene"] = body.source_gene
    _jobs[jid]["source_organism"] = body.source_organism
    _jobs[jid]["target_organism"] = body.target_organism
    _jobs[jid]["target_cell"] = body.target_cell

    async def _run_compat():
        try:
            from core.query_engine import run_compatibility_query  # noqa: WPS433
            result = await run_compatibility_query(
                source_gene=body.source_gene,
                source_organism=body.source_organism,
                target_organism=body.target_organism,
                target_cell=body.target_cell,
            )
            _jobs[jid].update({
                "status": "done",
                "compatibility_result": result.model_dump() if hasattr(result, "model_dump") else result,
            })
        except ImportError:
            _jobs[jid]["status"] = "error"
            _jobs[jid]["error"] = (
                "run_compatibility_query is not yet implemented in core.query_engine"
            )
        except Exception as exc:
            _jobs[jid]["status"] = "error"
            _jobs[jid]["error"] = str(exc)

    background_tasks.add_task(_run_compat)
    return JobCreatedResponse(job_id=jid)


@router.get("/jobs/{job_id}", response_model=JobStatusResponse, tags=["jobs"])
async def job_status(job_id: str):
    """Poll the status of any job (analyze or compare)."""
    from web.app import _jobs  # noqa: WPS433

    job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    return JobStatusResponse(
        job_id=job_id,
        status=job.get("status", "unknown"),
        created_at=job.get("created_at"),
        error=job.get("error"),
    )


@router.get("/jobs/{job_id}/report", tags=["jobs"])
async def job_report(job_id: str):
    """Return the full JSON payload for a completed job.

    Raises 404 if the job does not exist; 425 (Too Early) if still running.
    """
    from web.app import _jobs  # noqa: WPS433

    job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    status = job.get("status")
    if status in ("pending", "running"):
        raise HTTPException(status_code=425, detail="Job is still running")
    if status == "error":
        raise HTTPException(status_code=500, detail=job.get("error", "Unknown error"))

    return job
