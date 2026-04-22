from __future__ import annotations

from typing import Any

import requests


_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

_SEARCH_QUERY = """
query SearchTarget($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits {
      id
      name
      entity
      object {
        ... on Target {
          id
          approvedSymbol
          approvedName
          biotype
        }
      }
    }
  }
}
"""

_TARGET_QUERY = """
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    pathways {
      pathway
      pathwayId
      topLevelTerm
    }
    associatedDiseases(page: { index: 0, size: 20 }) {
      count
      rows {
        score
        disease {
          id
          name
          therapeuticAreas {
            name
          }
        }
        datasourceScores {
          id
          score
        }
      }
    }
    knownDrugs(size: 10) {
      count
      rows {
        drug {
          name
          maximumClinicalTrialPhase
          drugType
        }
        phase
        status
        disease {
          name
        }
        mechanismOfAction
      }
    }
  }
}
"""


class OpenTargetsLookup:
    """Fetch disease associations and drug data from the Open Targets Platform.

    Uses the public GraphQL API — no authentication required.
    """

    def __init__(self, *, timeout_seconds: float = 25.0, user_agent: str = "") -> None:
        self.timeout = timeout_seconds
        self.session = requests.Session()
        self.session.headers["Content-Type"] = "application/json"
        if user_agent:
            self.session.headers["User-Agent"] = user_agent

    def collect(self, *, gene_symbol: str) -> dict[str, Any]:
        if not gene_symbol:
            return {"open_targets_panel": None}

        ensembl_id = self._search_target(gene_symbol)
        if not ensembl_id:
            return {"open_targets_panel": None}

        data = self._fetch_target(ensembl_id)
        if not data:
            return {"open_targets_panel": None}

        target = data.get("target") or {}
        if not target:
            return {"open_targets_panel": None}

        assoc = target.get("associatedDiseases") or {}
        drugs_raw = target.get("knownDrugs") or {}

        disease_rows = self._disease_rows(assoc.get("rows", []))
        drug_rows = self._drug_rows(drugs_raw.get("rows", []))
        pathways = _pathways(target.get("pathways") or [])
        function_descriptions = (target.get("functionDescriptions") or [])[:3]

        panel: dict[str, Any] = {
            "ensembl_id": ensembl_id,
            "approved_name": target.get("approvedName", ""),
            "biotype": target.get("biotype", ""),
            "ot_url": f"https://platform.opentargets.org/target/{ensembl_id}",
            "total_disease_count": assoc.get("count", 0),
            "disease_rows": disease_rows,
            "total_drug_count": drugs_raw.get("count", 0),
            "drug_rows": drug_rows,
            "pathways": pathways,
            "function_descriptions": function_descriptions,
        }

        if not (disease_rows or drug_rows):
            return {"open_targets_panel": None}

        return {"open_targets_panel": panel}

    # -------------------------------------------------------------- #

    def _graphql(self, query: str, variables: dict) -> dict | None:
        try:
            resp = self.session.post(
                _GRAPHQL_URL,
                json={"query": query, "variables": variables},
                timeout=self.timeout,
            )
            resp.raise_for_status()
            result = resp.json()
            return result.get("data")
        except Exception:
            return None

    def _search_target(self, gene_symbol: str) -> str | None:
        data = self._graphql(_SEARCH_QUERY, {"q": gene_symbol})
        if not data:
            return None
        hits = (data.get("search") or {}).get("hits", [])
        for hit in hits:
            obj = hit.get("object") or {}
            symbol = obj.get("approvedSymbol", "")
            if symbol.upper() == gene_symbol.upper():
                return obj.get("id") or hit.get("id")
        # Fall back to first hit
        if hits:
            return (hits[0].get("object") or {}).get("id") or hits[0].get("id")
        return None

    def _fetch_target(self, ensembl_id: str) -> dict | None:
        return self._graphql(_TARGET_QUERY, {"ensemblId": ensembl_id})

    @staticmethod
    def _disease_rows(rows: list) -> list[dict[str, Any]]:
        result = []
        for row in rows:
            disease = row.get("disease") or {}
            score = row.get("score", 0.0)
            areas = [a["name"] for a in (disease.get("therapeuticAreas") or [])[:2]]
            result.append({
                "disease_name": disease.get("name", "—"),
                "disease_id": disease.get("id", ""),
                "score": round(score, 3),
                "score_pct": round(score * 100, 1),
                "therapeutic_areas": ", ".join(areas) if areas else "—",
                "ot_disease_url": f"https://platform.opentargets.org/disease/{disease.get('id', '')}",
            })
        result.sort(key=lambda r: r["score"], reverse=True)
        return result[:20]

    @staticmethod
    def _drug_rows(rows: list) -> list[dict[str, Any]]:
        result = []
        seen: set[str] = set()
        for row in rows:
            drug = row.get("drug") or {}
            name = drug.get("name", "")
            if not name or name in seen:
                continue
            seen.add(name)
            disease = row.get("disease") or {}
            result.append({
                "drug_name": name,
                "drug_type": drug.get("drugType", "—"),
                "max_phase": drug.get("maximumClinicalTrialPhase"),
                "phase_label": _phase_label(drug.get("maximumClinicalTrialPhase")),
                "mechanism": row.get("mechanismOfAction", "—"),
                "indication": disease.get("name", "—"),
            })
        return result


def _phase_label(phase: int | None) -> str:
    if phase is None:
        return "—"
    return {4: "Approved", 3: "Phase III", 2: "Phase II", 1: "Phase I"}.get(phase, f"Phase {phase}")


def _pathways(raw: list) -> list[dict[str, str]]:
    seen: set[str] = set()
    result = []
    for p in raw:
        name = p.get("pathway", "")
        pid = p.get("pathwayId", "")
        top = p.get("topLevelTerm", "")
        if name and name not in seen:
            seen.add(name)
            result.append({"name": name, "id": pid, "top_level": top})
    return result[:15]
