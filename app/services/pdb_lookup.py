from __future__ import annotations

from typing import Any

import requests


class PdbLookup:
    def __init__(
        self, *, base_url: str, timeout_seconds: float, user_agent: str = ""
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout_seconds
        self.session = requests.Session()
        if user_agent:
            self.session.headers["User-Agent"] = user_agent

    def collect(
        self, *, accession: str, variant_position: int
    ) -> dict[str, Any]:
        base_accession = accession.split("-")[0]
        entries = self._fetch_entries(base_accession, variant_position)
        if not entries:
            return {"pdb_panel": None}

        covering = [e for e in entries if e["covers_variant"]]
        best: dict[str, Any] | None = None
        experimental = [e for e in covering if e["method"] in {"X-ray", "EM", "Neutron"}]
        if experimental:
            best = min(
                experimental,
                key=lambda e: e["resolution_value"] if e["resolution_value"] is not None else 99.0,
            )
        elif covering:
            best = covering[0]

        return {
            "pdb_panel": {
                "entries": entries,
                "total_count": len(entries),
                "variant_covered_count": len(covering),
                "best_entry": best,
            }
        }

    def _fetch_entries(
        self, accession: str, variant_position: int
    ) -> list[dict[str, Any]]:
        try:
            resp = self.session.get(
                f"{self.base_url}/uniprotkb/{accession}.json",
                timeout=self.timeout,
            )
            resp.raise_for_status()
            payload = resp.json()
        except Exception:
            return []

        entries: list[dict[str, Any]] = []
        for ref in payload.get("uniProtKBCrossReferences", []):
            if ref.get("database") != "PDB":
                continue
            pdb_id = ref.get("id", "")
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            method = props.get("Method", "")
            resolution_str = props.get("Resolution", "")
            chains_str = props.get("Chains", "")

            resolution_value = _parse_resolution(resolution_str)
            chain_id, res_start, res_end = _parse_chains(chains_str)
            covers = (
                res_start is not None
                and res_end is not None
                and res_start <= variant_position <= res_end
            )

            entries.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "method": method,
                    "resolution": resolution_str,
                    "resolution_value": resolution_value,
                    "residue_start": res_start,
                    "residue_end": res_end,
                    "covers_variant": covers,
                    "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}",
                    "viewer_url": f"https://www.rcsb.org/3d-view/{pdb_id}",
                }
            )

        entries.sort(
            key=lambda e: (
                not e["covers_variant"],
                e["resolution_value"] if e["resolution_value"] is not None else 99.0,
            )
        )
        return entries


def _parse_resolution(res_str: str) -> float | None:
    if not res_str:
        return None
    try:
        return float(res_str.replace("A", "").strip().split()[0])
    except (ValueError, IndexError):
        return None


def _parse_chains(chains_str: str) -> tuple[str, int | None, int | None]:
    if not chains_str or "=" not in chains_str:
        return (chains_str, None, None)
    try:
        chain_part, range_part = chains_str.split("=", 1)
        start_str, end_str = range_part.split("-", 1)
        return (chain_part.strip(), int(start_str.strip()), int(end_str.strip()))
    except (ValueError, AttributeError):
        return (chains_str, None, None)
