from __future__ import annotations

import csv
import gzip
import io
import math
from typing import Any

import requests


PLDDT_LEVELS = [
    (90.0, "Very high", "plddt-very-high", "#0053D6"),
    (70.0, "High",      "plddt-high",      "#65CBF3"),
    (50.0, "Low",       "plddt-low",       "#FFDB13"),
    (0.0,  "Very low",  "plddt-very-low",  "#FF7D45"),
]

AM_CLASSES = {
    "likely_pathogenic": ("Likely pathogenic", "am-pathogenic"),
    "ambiguous":         ("Ambiguous",         "am-ambiguous"),
    "likely_benign":     ("Likely benign",     "am-benign"),
}

_API_BASE = "https://alphafold.ebi.ac.uk/api"
_VIEWER_BASE = "https://alphafold.ebi.ac.uk/entry"
_DOWNLOAD_LIMIT = 5 * 1024 * 1024  # 5 MB cap for PDB / AM files


class AlphaFoldLookup:
    def __init__(self, *, timeout_seconds: float = 20.0, user_agent: str = "") -> None:
        self.timeout = timeout_seconds
        self.session = requests.Session()
        if user_agent:
            self.session.headers["User-Agent"] = user_agent

    def collect(
        self,
        *,
        accession: str,
        variant_position: int,
        wild_type: str,
        mutant: str,
        sequence_length: int,
    ) -> dict[str, Any]:
        base_accession = accession.split("-")[0]
        meta = self._fetch_meta(base_accession)
        if meta is None:
            return {"alphafold_panel": None}

        residues = self._fetch_plddt(meta.get("pdbUrl", ""))
        am_score = self._fetch_am(
            meta.get("amAnnotationsUrl", ""),
            wild_type=wild_type,
            position=variant_position,
            mutant=mutant,
        )

        panel: dict[str, Any] = {
            "entry_id": meta.get("entryId", ""),
            "viewer_url": f"{_VIEWER_BASE}/{base_accession}",
            "model_date": meta.get("modelCreatedDate", ""),
            "uniprot_start": meta.get("uniprotStart", 1),
            "uniprot_end": meta.get("uniprotEnd", sequence_length),
            "plddt_available": bool(residues),
            "alphamissense": am_score,
        }

        if not residues:
            return {"alphafold_panel": panel}

        mean_plddt = sum(r["plddt"] for r in residues) / len(residues)
        counts: dict[str, int] = {}
        for r in residues:
            css = _plddt_css(r["plddt"])
            counts[css] = counts.get(css, 0) + 1

        variant_row = next(
            (
                {**r, **_plddt_meta(r["plddt"]), "height_pct": round(r["plddt"], 1)}
                for r in residues
                if r["position"] == variant_position
            ),
            None,
        )

        panel.update(
            {
                "mean_plddt": round(mean_plddt, 1),
                "mean_plddt_label": _plddt_label(mean_plddt),
                "mean_plddt_class": _plddt_css(mean_plddt),
                "variant_row": variant_row,
                "segments": _build_segments(residues, sequence_length),
                "confidence_cards": [
                    {
                        "label": label,
                        "class_name": css,
                        "color": color,
                        "count": counts.get(css, 0),
                        "percent": round((counts.get(css, 0) / len(residues)) * 100.0, 1),
                    }
                    for _, label, css, color in PLDDT_LEVELS
                ],
            }
        )
        return {"alphafold_panel": panel}

    # ------------------------------------------------------------------ #
    #  Private fetch helpers                                               #
    # ------------------------------------------------------------------ #

    def _fetch_meta(self, accession: str) -> dict[str, Any] | None:
        try:
            resp = self.session.get(
                f"{_API_BASE}/prediction/{accession}", timeout=self.timeout
            )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            data = resp.json()
            return data[0] if data else None
        except Exception:
            return None

    def _fetch_plddt(self, pdb_url: str) -> list[dict[str, Any]]:
        if not pdb_url:
            return []
        try:
            resp = self.session.get(
                pdb_url, timeout=max(self.timeout, 45.0), stream=True
            )
            resp.raise_for_status()
            chunks: list[bytes] = []
            total = 0
            for chunk in resp.iter_content(65536):
                chunks.append(chunk)
                total += len(chunk)
                if total >= _DOWNLOAD_LIMIT:
                    break
            text = b"".join(chunks).decode("utf-8", errors="replace")
            return _parse_pdb_plddt(text)
        except Exception:
            return []

    def _fetch_am(
        self, am_url: str, *, wild_type: str, position: int, mutant: str
    ) -> dict[str, Any] | None:
        if not am_url:
            return None
        target = f"{wild_type}{position}{mutant}"
        try:
            resp = self.session.get(
                am_url, timeout=max(self.timeout, 60.0), stream=True
            )
            resp.raise_for_status()
            chunks: list[bytes] = []
            total = 0
            for chunk in resp.iter_content(65536):
                chunks.append(chunk)
                total += len(chunk)
                if total >= 4 * 1024 * 1024:
                    break
            raw = b"".join(chunks)
            if am_url.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                except Exception:
                    pass
            content = raw.decode("utf-8", errors="replace")
            reader = csv.DictReader(io.StringIO(content))
            for row in reader:
                if row.get("protein_variant", "") == target:
                    score = float(row.get("am_pathogenicity", 0.0))
                    am_class = row.get("am_class", "")
                    label, css = AM_CLASSES.get(
                        am_class,
                        (am_class.replace("_", " ").title(), "am-unknown"),
                    )
                    return {
                        "score": round(score, 4),
                        "class": am_class,
                        "label": label,
                        "class_name": css,
                    }
        except Exception:
            pass
        return None


# ------------------------------------------------------------------ #
#  PDB parsing                                                         #
# ------------------------------------------------------------------ #


def _parse_pdb_plddt(pdb_text: str) -> list[dict[str, Any]]:
    residues: list[dict[str, Any]] = []
    seen: set[int] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            pos = int(line[22:26])
            plddt = float(line[60:66])
        except (ValueError, IndexError):
            continue
        if pos in seen:
            continue
        seen.add(pos)
        residues.append({"position": pos, "plddt": plddt})
    return sorted(residues, key=lambda r: r["position"])


# ------------------------------------------------------------------ #
#  pLDDT classification helpers                                        #
# ------------------------------------------------------------------ #


def _plddt_threshold(score: float) -> tuple[float, str, str, str]:
    for threshold, label, css, color in PLDDT_LEVELS:
        if score >= threshold:
            return threshold, label, css, color
    return 0.0, "Very low", "plddt-very-low", "#FF7D45"


def _plddt_label(score: float) -> str:
    return _plddt_threshold(score)[1]


def _plddt_css(score: float) -> str:
    return _plddt_threshold(score)[2]


def _plddt_meta(score: float) -> dict[str, str]:
    _, label, css, color = _plddt_threshold(score)
    return {"plddt_label": label, "plddt_class": css, "plddt_color": color}


# ------------------------------------------------------------------ #
#  Segment builder for the color strip                                 #
# ------------------------------------------------------------------ #


def _build_segments(
    residues: list[dict[str, Any]], sequence_length: int
) -> list[dict[str, Any]]:
    if not residues:
        return []
    segments: list[dict[str, Any]] = []
    current_css = _plddt_css(residues[0]["plddt"])
    seg_start = residues[0]["position"]
    seg_end = residues[0]["position"]

    for r in residues[1:]:
        css = _plddt_css(r["plddt"])
        if css == current_css:
            seg_end = r["position"]
        else:
            segments.append(
                _seg_payload(seg_start, seg_end, current_css, sequence_length)
            )
            current_css = css
            seg_start = seg_end = r["position"]

    segments.append(
        _seg_payload(seg_start, seg_end, current_css, sequence_length)
    )
    return segments


def _seg_payload(
    start: int, end: int, css: str, total: int
) -> dict[str, Any]:
    label = color = ""
    for _, lbl, c, clr in PLDDT_LEVELS:
        if c == css:
            label, color = lbl, clr
            break
    return {
        "start": start,
        "end": end,
        "class_name": css,
        "label": label,
        "color": color,
        "left_pct": round(((start - 1) / total) * 100.0, 3),
        "width_pct": max(round(((end - start + 1) / total) * 100.0, 3), 0.5),
    }
