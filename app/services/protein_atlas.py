from __future__ import annotations

from typing import Any

import requests


_API_BASE = "https://www.proteinatlas.org"


class ProteinAtlasLookup:
    """Fetch expression, subcellular localisation, and disease data
    from the Human Protein Atlas public JSON search API.

    Endpoint: /search/{gene_symbol}?format=json  (returns a JSON array)
    Each item is a flat dict of summary fields — no per-tissue rows.
    """

    def __init__(self, *, timeout_seconds: float = 20.0, user_agent: str = "") -> None:
        self.timeout = timeout_seconds
        self.session = requests.Session()
        if user_agent:
            self.session.headers["User-Agent"] = user_agent

    def collect(self, *, gene_symbol: str) -> dict[str, Any]:
        if not gene_symbol:
            return {"protein_atlas_panel": None}

        data = self._fetch(gene_symbol)
        if data is None:
            return {"protein_atlas_panel": None}

        subcell = self._subcellular(data)
        diseases = self._diseases(data)
        bio_proc = self._list_field(data, "Biological process")
        mol_func = self._list_field(data, "Molecular function")
        protein_class = self._list_field(data, "Protein class")
        rna_fields = self._rna_summary(data)

        panel: dict[str, Any] = {
            "gene_symbol": gene_symbol,
            "hpa_url": f"{_API_BASE}/{gene_symbol}",
            "ensembl_id": data.get("Ensembl", ""),
            "gene_description": data.get("Gene description", ""),
            "protein_class": protein_class,
            "biological_process": bio_proc,
            "molecular_function": mol_func,
            "subcellular_locations": subcell,
            "disease_associations": diseases,
            "rna_fields": rna_fields,
            "interactions": data.get("Interactions", ""),
            "evidence": data.get("Evidence", ""),
            "tissue_expression_cluster": data.get("Tissue expression cluster", ""),
        }

        if not (subcell or diseases or protein_class or rna_fields):
            return {"protein_atlas_panel": None}

        return {"protein_atlas_panel": panel}

    # -------------------------------------------------------------- #

    def _fetch(self, gene_symbol: str) -> dict[str, Any] | None:
        try:
            resp = self.session.get(
                f"{_API_BASE}/search/{gene_symbol}",
                params={"format": "json"},
                timeout=self.timeout,
            )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            result = resp.json()
            if isinstance(result, list) and result:
                for item in result:
                    gene = (item.get("Gene") or "").upper()
                    if gene == gene_symbol.upper():
                        return item
                return result[0]
            if isinstance(result, dict):
                return result
            return None
        except Exception:
            return None

    @staticmethod
    def _list_field(data: dict, key: str) -> list[str]:
        val = data.get(key, [])
        if isinstance(val, list):
            return [str(v).strip() for v in val if v and str(v).strip() not in ("None", "")]
        if isinstance(val, str) and val and val != "None":
            return [s.strip() for s in val.split(",") if s.strip()]
        return []

    @staticmethod
    def _subcellular(data: dict) -> list[dict[str, str]]:
        """Subcellular locations from main + additional location fields."""
        locs: list[str] = []
        for key in ("Subcellular main location", "Subcellular additional location", "Subcellular location"):
            val = data.get(key, [])
            if isinstance(val, list):
                locs.extend([str(v).strip() for v in val if v and str(v).strip() not in ("None", "")])
            elif isinstance(val, str) and val and val != "None":
                locs.extend([s.strip() for s in val.split(",") if s.strip()])
        seen: set[str] = set()
        result = []
        for loc in locs:
            if loc and loc not in seen:
                seen.add(loc)
                result.append({"location": loc, "reliability": ""})
        return result

    @staticmethod
    def _diseases(data: dict) -> list[dict[str, str]]:
        raw = data.get("Disease involvement", [])
        if isinstance(raw, list):
            return [{"disease": str(d).strip(), "source": "HPA"} for d in raw
                    if d and str(d).strip() not in ("None", "")]
        if isinstance(raw, str) and raw and raw != "None":
            return [{"disease": d.strip(), "source": "HPA"} for d in raw.split(",") if d.strip()]
        return []

    @staticmethod
    def _rna_summary(data: dict) -> list[dict[str, str]]:
        """Pull out key RNA expression summary strings."""
        fields = [
            ("Tissue specificity", "RNA tissue specificity"),
            ("Tissue distribution", "RNA tissue distribution"),
            ("Single cell specificity", "RNA single cell type specificity"),
            ("Single cell distribution", "RNA single cell type distribution"),
            ("Blood cell specificity", "RNA blood cell specificity"),
            ("Blood cell distribution", "RNA blood cell distribution"),
            ("Cancer specificity", "RNA cancer specificity"),
            ("Cancer distribution", "RNA cancer distribution"),
        ]
        result = []
        for label, key in fields:
            val = data.get(key, "")
            if val and str(val) not in ("None", "null", ""):
                result.append({"label": label, "value": str(val)})
        return result
