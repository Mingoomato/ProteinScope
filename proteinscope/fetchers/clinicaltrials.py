"""ClinicalTrials.gov v2 REST API fetcher.

Source: Clinicaltrials.gov API v2 (Vivli, 2023): https://clinicaltrials.gov/data-api/api
No API key required. Returns empty list on any failure.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

_FIELDS = "nctId,briefTitle,phase,overallStatus,interventionType,sponsorName"


async def fetch_clinical_trials(gene_name: str) -> list[dict]:
    """Fetch ClinicalTrials.gov v2 studies mentioning *gene_name*.

    Args:
        gene_name: HGNC gene symbol to search (e.g. "EGFR", "BRCA1").

    Returns:
        List of dicts with keys:
            nct_id, title, phase, status, intervention_type, sponsor
        Empty list on any network or parse failure — never propagates.

    Source:
        ClinicalTrials.gov REST API v2 (Vivli, 2023).
        https://clinicaltrials.gov/data-api/api
    """
    import httpx

    url = (
        f"{_BASE_URL}"
        f"?query.term={gene_name}"
        f"&fields={_FIELDS}"
        f"&pageSize=20"
    )

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            response = await client.get(url)
            response.raise_for_status()
            payload = response.json()

        studies: list[dict] = payload.get("studies", [])
        results: list[dict] = []

        for study in studies:
            protocol = study.get("protocolSection", {})

            # identificationModule
            id_module = protocol.get("identificationModule", {})
            nct_id: str = id_module.get("nctId", "")
            title: str = id_module.get("briefTitle", "")

            # designModule → phases list (e.g. ["PHASE2", "PHASE3"])
            design_module = protocol.get("designModule", {})
            phases: list[str] = design_module.get("phases", [])
            phase: str = ", ".join(phases) if phases else "N/A"

            # statusModule
            status_module = protocol.get("statusModule", {})
            overall_status: str = status_module.get("overallStatus", "")

            # armsInterventionsModule → interventions list
            arms_module = protocol.get("armsInterventionsModule", {})
            interventions: list[dict] = arms_module.get("interventions", [])
            intervention_types: list[str] = list(
                {iv.get("type", "") for iv in interventions if iv.get("type")}
            )
            intervention_type: str = ", ".join(intervention_types) if intervention_types else "N/A"

            # sponsorCollaboratorsModule
            sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
            lead_sponsor: str = (
                sponsor_module.get("leadSponsor", {}).get("name", "")
            )

            if not nct_id:
                continue

            results.append(
                {
                    "nct_id": nct_id,
                    "title": title,
                    "phase": phase,
                    "status": overall_status,
                    "intervention_type": intervention_type,
                    "sponsor": lead_sponsor,
                }
            )

        return results

    except Exception as exc:  # noqa: BLE001
        logger.debug("ClinicalTrials fetch failed for %r: %s", gene_name, exc)
        return []
