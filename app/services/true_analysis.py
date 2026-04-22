from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import shlex
import subprocess
import uuid
from typing import Any

from app.services.tool_registry import PREDICTOR_SPECS, build_resource_cards, resolve_command_template
from app.services.variant_parser import MissenseVariant
from app.services.protein_lookup import ProteinRecord


CATEGORY_LABELS = {
    "domains": "Domains and motifs",
    "structure": "Structure and accessibility",
    "topology": "Topology and targeting",
    "localization": "Localization and sorting",
    "disorder": "Disorder and flexibility",
    "ptm": "PTM predictions",
    "conservation": "Conservation and coupling",
}

CATEGORY_COLORS = {
    "domains": "#2c8f71",
    "structure": "#4b6cb7",
    "topology": "#d86f45",
    "localization": "#c08a1f",
    "disorder": "#8d5bd1",
    "ptm": "#b53f66",
    "conservation": "#3d5f9d",
}


@dataclass(frozen=True)
class PredictorStatus:
    source: str
    status: str
    summary: str


class TrueAnalysisError(RuntimeError):
    """Raised when a true-analysis tool contract fails."""


class TrueAnalysisOrchestrator:
    def __init__(self, *, workdir: str, timeout_seconds: float, enabled_predictors: set[str] | None = None) -> None:
        self.workdir = Path(workdir)
        self.timeout_seconds = timeout_seconds
        self.enabled_predictors = enabled_predictors or {spec.key for spec in PREDICTOR_SPECS}
        self.workdir.mkdir(parents=True, exist_ok=True)

    def collect(self, protein: ProteinRecord, variant: MissenseVariant) -> dict[str, Any]:
        interval_features: list[dict[str, Any]] = []
        site_features: list[dict[str, Any]] = []
        evidence_rows: list[dict[str, Any]] = []
        residue_tracks: list[dict[str, Any]] = []
        predictor_tables: list[dict[str, Any]] = []
        predictor_statuses: list[dict[str, Any]] = []
        notices: list[str] = []

        resource_cards = build_resource_cards(self.enabled_predictors)

        for spec in PREDICTOR_SPECS:
            if spec.key not in self.enabled_predictors:
                continue
            command_template = resolve_command_template(spec)
            if not command_template:
                predictor_statuses.append(
                    PredictorStatus(
                        source=spec.label,
                        status="missing",
                        summary=f"Set {spec.command_env} to enable this predictor. {spec.resource_hint}",
                    ).__dict__
                )
                continue

            try:
                result = self._run_predictor(spec, command_template, protein, variant)
                interval_features.extend(result.get("interval_features", []))
                site_features.extend(result.get("site_features", []))
                evidence_rows.extend(result.get("evidence_rows", []))
                residue_tracks.extend(result.get("residue_tracks", []))
                predictor_tables.extend(result.get("tables", []))
                notices.extend(result.get("notices", []))
                predictor_statuses.append(
                    PredictorStatus(
                        source=spec.label,
                        status="ok",
                        summary=result.get("summary") or f"{spec.label} completed successfully.",
                    ).__dict__
                )
            except TrueAnalysisError as exc:
                predictor_statuses.append(
                    PredictorStatus(
                        source=spec.label,
                        status="error",
                        summary=str(exc),
                    ).__dict__
                )

        architecture_rows = self._build_architecture_rows(interval_features, protein.length, variant.position)
        feature_counts = self._build_feature_counts(interval_features, site_features)
        evidence_plot = self._build_evidence_plot(evidence_rows, protein.length, variant.position)

        return {
            "interval_features": sorted(interval_features, key=lambda row: (row["start"], row["end"], row["label"])),
            "site_features": sorted(site_features, key=lambda row: (row["position"], row["label"])),
            "evidence_rows": evidence_rows,
            "residue_tracks": residue_tracks,
            "predictor_tables": predictor_tables,
            "resource_cards": resource_cards,
            "architecture_rows": architecture_rows,
            "feature_counts": feature_counts,
            "evidence_plot": evidence_plot,
            "source_status": predictor_statuses,
            "notices": sorted(set(notices)),
        }

    def _run_predictor(
        self,
        spec,
        command_template: str,
        protein: ProteinRecord,
        variant: MissenseVariant,
    ) -> dict[str, Any]:
        run_dir = self.workdir / f"{spec.key}-{uuid.uuid4().hex[:12]}"
        output_dir = run_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        input_fasta = run_dir / "input.fasta"
        output_json = output_dir / "result.json"

        input_fasta.write_text(
            f">{protein.accession}|{protein.gene_symbol or protein.entry_name}|{variant.short_notation}\n{protein.sequence}\n",
            encoding="utf-8",
        )

        command = command_template.format(
            input_fasta=str(input_fasta),
            output_dir=str(output_dir),
            output_json=str(output_json),
            sequence_id=protein.accession,
            gene_symbol=protein.gene_symbol or protein.entry_name,
            variant=variant.short_notation,
            variant_position=variant.position,
            wild_type=variant.wild_type,
            mutant=variant.mutant,
            organism="human",
            workdir=str(run_dir),
        )

        try:
            completed = subprocess.run(
                shlex.split(command),
                cwd=run_dir,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds,
                check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise TrueAnalysisError(f"{spec.label} failed to start: {exc}") from exc

        if completed.returncode != 0:
            stderr = (completed.stderr or completed.stdout or "").strip()
            raise TrueAnalysisError(f"{spec.label} command failed: {stderr or 'unknown error'}")

        if not output_json.exists():
            raise TrueAnalysisError(
                f"{spec.label} completed, but the wrapper did not create {output_json.name}."
            )

        try:
            payload = json.loads(output_json.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise TrueAnalysisError(f"{spec.label} wrapper wrote invalid JSON output.") from exc

        normalized = self._normalize_payload(spec, payload, variant.position)
        if completed.stderr.strip():
            normalized.setdefault("notices", []).append(completed.stderr.strip())
        return normalized

    def _normalize_payload(self, spec, payload: dict[str, Any], variant_position: int) -> dict[str, Any]:
        normalized = {
            "summary": payload.get("summary", ""),
            "interval_features": [],
            "site_features": [],
            "evidence_rows": [],
            "residue_tracks": [],
            "tables": payload.get("tables", []),
            "notices": payload.get("notices", []),
        }

        for feature in payload.get("interval_features", []):
            normalized["interval_features"].append(
                {
                    "source": feature.get("source", spec.label),
                    "type": feature.get("type", spec.key.upper()),
                    "category": feature.get("category", spec.category),
                    "label": feature.get("label", spec.label),
                    "description": feature.get("description", feature.get("label", spec.label)),
                    "start": int(feature["start"]),
                    "end": int(feature["end"]),
                    "evidence": feature.get("evidence", ""),
                    "variant_hit": bool(feature.get("variant_hit", False)),
                }
            )

        for feature in payload.get("site_features", []):
            position = int(feature["position"])
            normalized["site_features"].append(
                {
                    "source": feature.get("source", spec.label),
                    "type": feature.get("type", spec.key.upper()),
                    "category": feature.get("category", spec.category),
                    "label": feature.get("label", spec.label),
                    "description": feature.get("description", feature.get("label", spec.label)),
                    "position": position,
                    "start": position,
                    "end": position,
                    "evidence": feature.get("evidence", ""),
                    "variant_hit": bool(feature.get("variant_hit", False)),
                }
            )

        for row in payload.get("evidence_rows", []):
            normalized["evidence_rows"].append(
                {
                    "source": row.get("source", spec.label),
                    "kind": row.get("kind", spec.key),
                    "label": row.get("label", spec.label),
                    "position": int(row["position"]),
                    "summary": row.get("summary", ""),
                    "link_hint": row.get("link_hint", ""),
                    "score": row.get("score", ""),
                }
            )

        for track in payload.get("residue_tracks", []):
            values = []
            for point in track.get("values", []):
                position = int(point["position"])
                raw_value = float(point["value"])
                values.append(
                    {
                        "position": position,
                        "value": raw_value,
                        "height_pct": max(min(raw_value * 100.0, 100.0), 0.0),
                        "is_variant_site": bool(point.get("is_variant_site", position == variant_position)),
                    }
                )
            normalized["residue_tracks"].append(
                {
                    "key": track.get("key", spec.key),
                    "label": track.get("label", spec.label),
                    "source": spec.label,
                    "values": values,
                }
            )

        return normalized

    @staticmethod
    def _build_architecture_rows(
        interval_features: list[dict[str, Any]], sequence_length: int, variant_position: int
    ) -> list[dict[str, Any]]:
        grouped: dict[str, list[dict[str, Any]]] = {}
        for feature in interval_features:
            grouped.setdefault(feature["category"], []).append(feature)

        rows = []
        for category in ("domains", "structure", "topology", "localization", "disorder", "ptm", "conservation"):
            items = grouped.get(category, [])
            segments = []
            for feature in items[:12]:
                segments.append(
                    {
                        **feature,
                        "left_pct": round(((feature["start"] - 1) / sequence_length) * 100, 3),
                        "width_pct": max(round(((feature["end"] - feature["start"] + 1) / sequence_length) * 100, 3), 0.7),
                        "color": CATEGORY_COLORS[category],
                    }
                )
            if not segments:
                continue
            rows.append(
                {
                    "key": category,
                    "label": CATEGORY_LABELS[category],
                    "segments": segments,
                    "variant_pct": round((variant_position / sequence_length) * 100, 3),
                }
            )
        return rows

    @staticmethod
    def _build_feature_counts(
        interval_features: list[dict[str, Any]], site_features: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        counts: dict[str, int] = {}
        for feature in interval_features + site_features:
            label = CATEGORY_LABELS.get(feature["category"], feature["category"].title())
            counts[label] = counts.get(label, 0) + 1
        return [{"label": label, "count": count} for label, count in sorted(counts.items())]

    @staticmethod
    def _build_evidence_plot(
        evidence_rows: list[dict[str, Any]], sequence_length: int, variant_position: int
    ) -> list[dict[str, Any]]:
        points = []
        for row in evidence_rows:
            points.append(
                {
                    **row,
                    "left_pct": round((row["position"] / sequence_length) * 100, 3),
                    "hit_variant_site": row["position"] == variant_position,
                }
            )
        return points
