"""OpenTargets Platform API fetcher — disease association scores.

Uses the GraphQL API (v4, free, no key required).
Requires an Ensembl gene ID (ENSGXXXXXXXXXXX) obtained from UniProt dbReferences.

# Citation: 김박사 (Codex GPT-5.4), Grand Consortium V3 2026-03-20
# API docs: https://platform-docs.opentargets.org/data-access/graphql-api
"""

from __future__ import annotations

import httpx

_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# GraphQL query confirmed by 김박사 via OpenTargets community docs:
# https://community.opentargets.org/t/returning-all-associations-data-using-the-platform-api/324
_QUERY = """
query targetDiseaseAssociations($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
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
"""


async def fetch_opentargets(ensembl_id: str) -> list[dict]:
    """Return ranked disease associations from OpenTargets for an Ensembl gene ID.

    Args:
        ensembl_id: Ensembl gene ID, e.g. "ENSG00000146648" (EGFR).
                    Obtained from UniProt entry["dbReferences"] where type=="Ensembl".

    Returns:
        List of dicts: [{disease_id, disease_name, score}]
        Empty list on any error or missing ID.
    """
    if not ensembl_id or not ensembl_id.startswith("ENSG"):
        return []

    try:
        async with httpx.AsyncClient(timeout=30) as client:
            resp = await client.post(
                _GRAPHQL_URL,
                json={"query": _QUERY, "variables": {"ensemblId": ensembl_id}},
                headers={"Content-Type": "application/json"},
            )
        if resp.status_code != 200:
            return []

        data = resp.json()
        target = (data.get("data") or {}).get("target") or {}
        rows = (target.get("associatedDiseases") or {}).get("rows") or []

        results = []
        for row in rows:
            disease = row.get("disease") or {}
            disease_id = disease.get("id", "")
            disease_name = disease.get("name", "")
            score = float(row.get("score") or 0.0)
            if disease_id and disease_name:
                results.append({
                    "disease_id": disease_id,
                    "disease_name": disease_name,
                    "score": score,
                    "source": "OpenTargets",
                })
        return results

    except Exception:
        return []


def extract_ensembl_id_from_uniprot(entry: dict) -> str | None:
    """Extract the first Ensembl gene ID from a UniProt entry dict.

    UniProt entry["dbReferences"] contains cross-references to Ensembl.
    Returns the gene-level ID (ENSG...), not transcript (ENST...).
    """
    for ref in entry.get("dbReferences", []):
        if ref.get("type") == "Ensembl":
            ref_id = ref.get("id", "")
            if ref_id.startswith("ENSG"):
                return ref_id
            # Some Ensembl refs have gene ID in properties
            for prop in ref.get("properties", []):
                if prop.get("key") == "gene ID":
                    gene_id = prop.get("value", "")
                    if gene_id.startswith("ENSG"):
                        return gene_id
    return None
