from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isfinite, sin, cos, sqrt, radians
from typing import Any

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:  # pragma: no cover - runtime dependency check
    ProteinAnalysis = None

from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import MissenseVariant


HYDROPATHY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

RESIDUE_MASSES = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.15,
    "Q": 146.15,
    "E": 147.13,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "L": 131.17,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

POLAR = set("RNDQEKHSTY")
POSITIVE = set("RKH")
NEGATIVE = set("DE")
AROMATIC = set("FWYH")

PKA = {
    "Cterm": 3.1,
    "Nterm": 8.0,
    "C": 8.5,
    "D": 3.9,
    "E": 4.1,
    "H": 6.5,
    "K": 10.8,
    "R": 12.5,
    "Y": 10.1,
}

STRUCTURE_STATE_META = {
    "H": {"label": "Alpha helix", "class_name": "helix", "accent": "#c65d44"},
    "E": {"label": "Beta strand", "class_name": "strand", "accent": "#356d9c"},
    "C": {"label": "Coil", "class_name": "coil", "accent": "#6d7f69"},
}


class AnalysisError(RuntimeError):
    """Raised when mutation analysis cannot be completed."""


@dataclass(frozen=True)
class MetricDelta:
    label: str
    wild_type: float | None
    mutant: float | None
    delta: float | None


class SequenceAnalysisService:
    def __init__(self, analysis_sections: list[dict[str, Any]]) -> None:
        self.analysis_sections = analysis_sections

    def build_report(
        self,
        *,
        protein: ProteinRecord,
        variant: MissenseVariant,
        annotations: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        self._validate_variant(protein.sequence, variant)
        mutant_sequence = self._mutate_sequence(protein.sequence, variant)
        metrics = self._metric_deltas(protein.sequence, mutant_sequence)
        context = self._residue_context(protein.sequence, variant)

        annotation_bundle = annotations or {
            "interval_features": [],
            "site_features": [],
            "evidence_rows": [],
            "residue_tracks": [],
            "predictor_tables": [],
            "resource_cards": [],
            "architecture_rows": [],
            "feature_counts": [],
            "evidence_plot": [],
            "source_status": [],
            "notices": [],
        }

        predictor_tables = annotation_bundle.get("predictor_tables", [])
        residue_tracks = annotation_bundle.get("residue_tracks", [])

        overlapping_intervals = [
            feature
            for feature in annotation_bundle["interval_features"]
            if feature["start"] <= variant.position <= feature["end"]
        ]
        overlapping_sites = [
            feature for feature in annotation_bundle["site_features"] if feature["position"] == variant.position
        ]
        active_source_status = self._active_source_status(annotation_bundle.get("source_status", []))
        active_resource_cards = self._active_resource_cards(
            annotation_bundle.get("resource_cards", []), active_source_status
        )
        accessibility_panel = self._accessibility_panel(predictor_tables, variant.position)
        secondary_structure_panel = self._secondary_structure_panel(
            predictor_tables, variant.position, protein.length
        )
        structure_interval_features = [
            feature for feature in annotation_bundle.get("interval_features", []) if feature.get("category") == "structure"
        ]
        structure_site_features = [
            feature for feature in annotation_bundle.get("site_features", []) if feature.get("category") == "structure"
        ]
        structure_architecture_rows = [
            row for row in annotation_bundle.get("architecture_rows", []) if row.get("key") == "structure"
        ]

        structural_verdict = self._structural_verdict(
            secondary_structure_panel=secondary_structure_panel,
            accessibility_panel=accessibility_panel,
            overlapping_intervals=overlapping_intervals,
        )
        amphipathic_panel = self._amphipathic_helix_panel(
            secondary_structure_panel=secondary_structure_panel,
            sequence=protein.sequence,
            variant_position=variant.position,
        )

        report = {
            "protein": protein,
            "variant": variant,
            "mutant_sequence": mutant_sequence,
            "context": context,
            "metric_deltas": metrics,
            "metric_plot": self._metric_plot(metrics),
            "hydropathy_window_plot": self._hydropathy_window_plot(context),
            "section_status": self.analysis_sections,
            "interval_features": annotation_bundle.get("interval_features", []),
            "site_features": annotation_bundle.get("site_features", []),
            "evidence_rows": annotation_bundle.get("evidence_rows", []),
            "residue_tracks": residue_tracks,
            "predictor_tables": predictor_tables,
            "supplementary_residue_tracks": self._supplementary_residue_tracks(residue_tracks),
            "supplementary_predictor_tables": self._supplementary_predictor_tables(predictor_tables),
            "resource_cards": active_resource_cards,
            "architecture_rows": annotation_bundle.get("architecture_rows", []),
            "structure_architecture_rows": structure_architecture_rows,
            "feature_counts": annotation_bundle.get("feature_counts", []),
            "evidence_plot": annotation_bundle.get("evidence_plot", []),
            "source_status": active_source_status,
            "notices": annotation_bundle.get("notices", []),
            "source_provenance": self._source_provenance(active_source_status),
            "summary_cards": self._summary_cards(protein, variant, metrics, overlapping_intervals, overlapping_sites),
            "accessibility_panel": accessibility_panel,
            "secondary_structure_panel": secondary_structure_panel,
            "structure_interval_features": structure_interval_features,
            "structure_site_features": structure_site_features,
            "structural_verdict": structural_verdict,
            "amphipathic_panel": amphipathic_panel,
        }
        return report

    @staticmethod
    def _validate_variant(sequence: str, variant: MissenseVariant) -> None:
        if variant.position > len(sequence):
            raise AnalysisError(
                f"Residue {variant.position} is outside the resolved sequence length of {len(sequence)}."
            )

        observed = sequence[variant.position - 1]
        if observed != variant.wild_type:
            raise AnalysisError(
                f"Wild-type residue mismatch at position {variant.position}: sequence has {observed}, "
                f"but the submitted variant expects {variant.wild_type}."
            )

    @staticmethod
    def _mutate_sequence(sequence: str, variant: MissenseVariant) -> str:
        index = variant.position - 1
        return sequence[:index] + variant.mutant + sequence[index + 1 :]

    def _metric_deltas(self, wild_type_sequence: str, mutant_sequence: str) -> list[MetricDelta]:
        wt_metrics = self._sequence_metrics(wild_type_sequence)
        mut_metrics = self._sequence_metrics(mutant_sequence)

        deltas = []
        for label in (
            "Molecular weight",
            "Isoelectric point",
            "GRAVY",
            "Instability index",
            "Aromaticity",
            "Aliphatic index",
        ):
            wt = wt_metrics.get(label)
            mut = mut_metrics.get(label)
            delta = round(mut - wt, 4) if wt is not None and mut is not None else None
            deltas.append(MetricDelta(label=label, wild_type=wt, mutant=mut, delta=delta))
        return deltas

    def _sequence_metrics(self, sequence: str) -> dict[str, float | None]:
        if ProteinAnalysis is not None:
            analysis = ProteinAnalysis(sequence)
            return {
                "Molecular weight": round(analysis.molecular_weight(), 3),
                "Isoelectric point": round(analysis.isoelectric_point(), 3),
                "GRAVY": round(analysis.gravy(), 3),
                "Instability index": round(analysis.instability_index(), 3),
                "Aromaticity": round(analysis.aromaticity(), 3),
                "Aliphatic index": round(self._aliphatic_index(sequence), 3),
            }

        return {
            "Molecular weight": round(self._approx_molecular_weight(sequence), 3),
            "Isoelectric point": round(self._approx_isoelectric_point(sequence), 3),
            "GRAVY": round(self._average_hydropathy(sequence), 3),
            "Instability index": None,
            "Aromaticity": round(self._aromaticity(sequence), 3),
            "Aliphatic index": round(self._aliphatic_index(sequence), 3),
        }

    @staticmethod
    def _approx_molecular_weight(sequence: str) -> float:
        total = sum(RESIDUE_MASSES[residue] for residue in sequence)
        water_loss = max(len(sequence) - 1, 0) * 18.015
        return total - water_loss

    @staticmethod
    def _approx_isoelectric_point(sequence: str) -> float:
        low, high = 0.0, 14.0
        for _ in range(60):
            mid = (low + high) / 2
            charge = _net_charge(sequence, mid)
            if charge > 0:
                low = mid
            else:
                high = mid
        return (low + high) / 2

    @staticmethod
    def _average_hydropathy(sequence: str) -> float:
        return sum(HYDROPATHY[residue] for residue in sequence) / len(sequence)

    @staticmethod
    def _aromaticity(sequence: str) -> float:
        aromatic_count = sum(1 for residue in sequence if residue in {"F", "W", "Y"})
        return aromatic_count / len(sequence)

    @staticmethod
    def _aliphatic_index(sequence: str) -> float:
        counts = Counter(sequence)
        length = len(sequence)
        if not length:
            return 0.0
        return ((counts["A"] + 2.9 * counts["V"] + 3.9 * (counts["I"] + counts["L"])) / length) * 100.0

    @staticmethod
    def _residue_context(sequence: str, variant: MissenseVariant) -> dict[str, Any]:
        window_radius = 10
        center_index = variant.position - 1
        start = max(0, center_index - window_radius)
        end = min(len(sequence), center_index + window_radius + 1)

        wild_window = sequence[start:end]
        mutant_window = wild_window[: center_index - start] + variant.mutant + wild_window[center_index - start + 1 :]
        reference = variant.wild_type
        alternate = variant.mutant

        return {
            "window_start": start + 1,
            "window_end": end,
            "wild_window": wild_window,
            "mutant_window": mutant_window,
            "window_variant_index": center_index - start,
            "wild_type_properties": {
                "hydropathy": HYDROPATHY[reference],
                "charge_class": _charge_class(reference),
                "polarity_class": _polarity_class(reference),
                "aromatic": reference in AROMATIC,
            },
            "mutant_properties": {
                "hydropathy": HYDROPATHY[alternate],
                "charge_class": _charge_class(alternate),
                "polarity_class": _polarity_class(alternate),
                "aromatic": alternate in AROMATIC,
            },
        }

    @staticmethod
    def _metric_plot(metrics: list[MetricDelta]) -> list[dict[str, Any]]:
        finite_deltas = [abs(metric.delta) for metric in metrics if metric.delta is not None and isfinite(metric.delta)]
        max_delta = max(finite_deltas, default=1.0)
        rows = []
        for metric in metrics:
            if metric.delta is None:
                rows.append(
                    {
                        "label": metric.label,
                        "width_pct": 0,
                        "direction": "missing",
                        "display_delta": "Unavailable",
                    }
                )
                continue
            rows.append(
                {
                    "label": metric.label,
                    "width_pct": round((abs(metric.delta) / max_delta) * 100, 2) if max_delta else 0,
                    "direction": "up" if metric.delta >= 0 else "down",
                    "display_delta": f"{metric.delta:+.3f}",
                }
            )
        return rows

    @staticmethod
    def _hydropathy_window_plot(context: dict[str, Any]) -> list[dict[str, Any]]:
        wild_window = context["wild_window"]
        mutant_window = context["mutant_window"]
        values = []
        scale = max(abs(score) for score in HYDROPATHY.values())
        for index, (wild, mutant) in enumerate(zip(wild_window, mutant_window), start=context["window_start"]):
            wild_score = HYDROPATHY[wild]
            mutant_score = HYDROPATHY[mutant]
            values.append(
                {
                    "position": index,
                    "wild_residue": wild,
                    "mutant_residue": mutant,
                    "wild_score": wild_score,
                    "mutant_score": mutant_score,
                    "wild_height_pct": round((abs(wild_score) / scale) * 100, 2),
                    "mutant_height_pct": round((abs(mutant_score) / scale) * 100, 2),
                    "wild_direction": "positive" if wild_score >= 0 else "negative",
                    "mutant_direction": "positive" if mutant_score >= 0 else "negative",
                    "is_variant_site": index == context["window_start"] + context["window_variant_index"],
                }
            )
        return values

    @staticmethod
    def _summary_cards(
        protein: ProteinRecord,
        variant: MissenseVariant,
        metrics: list[MetricDelta],
        overlapping_intervals: list[dict[str, Any]],
        overlapping_sites: list[dict[str, Any]],
    ) -> list[dict[str, str]]:
        hydropathy_metric = next(metric for metric in metrics if metric.label == "GRAVY")
        interval_text = (
            ", ".join(feature["label"] for feature in overlapping_intervals[:3]) if overlapping_intervals else "No interval feature overlap found"
        )
        site_text = (
            ", ".join(feature["label"] for feature in overlapping_sites[:3]) if overlapping_sites else "No exact site feature overlap found"
        )
        return [
            {
                "label": "Target scope",
                "value": f"{protein.gene_symbol or protein.entry_name} / {protein.accession}",
            },
            {
                "label": "Feature overlap",
                "value": interval_text,
            },
            {
                "label": "Residue-site overlap",
                "value": site_text,
            },
            {
                "label": "Hydropathy delta",
                "value": f"{hydropathy_metric.delta:+.3f}" if hydropathy_metric.delta is not None else "Unavailable",
            },
            {
                "label": "Neighborhood",
                "value": f"Residues {max(1, variant.position - 10)}-{min(protein.length, variant.position + 10)}",
            },
        ]

    @staticmethod
    def _source_provenance(source_status: list[dict[str, Any]]) -> list[dict[str, str]]:
        rows = [
            {
                "source": "UniProt REST",
                "purpose": "Human protein resolution and canonical sequence retrieval for the analysis input.",
            },
            {
                "source": "Internal sequence metrics",
                "purpose": "Wild-type versus mutant property deltas and residue-window plots.",
            },
        ]
        for status in source_status:
            rows.append(
                {
                    "source": status["source"],
                    "purpose": status["summary"],
                }
            )
        return rows

    @staticmethod
    def _active_source_status(source_status: list[dict[str, Any]]) -> list[dict[str, Any]]:
        return [status for status in source_status if status.get("status") == "ok"]

    @staticmethod
    def _active_resource_cards(
        resource_cards: list[dict[str, Any]], source_status: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        active_labels = {status.get("source") for status in source_status}
        return [card for card in resource_cards if card.get("label") in active_labels]

    @staticmethod
    def _accessibility_panel(predictor_tables: list[dict[str, Any]], variant_position: int) -> dict[str, Any] | None:
        residue_table = next(
            (table for table in predictor_tables if table.get("title") == "DMVFL-RSA residue accessibility"),
            None,
        )
        if residue_table is None:
            return None

        rows = residue_table.get("rows", [])
        if not rows:
            return None

        max_asa = max(float(row["ASA"]) for row in rows) or 1.0
        mean_rsa = sum(float(row["RSA"]) for row in rows) / len(rows)
        mean_asa = sum(float(row["ASA"]) for row in rows) / len(rows)
        exposed_count = sum(1 for row in rows if row["Exposure"] == "Exposed")
        buried_count = sum(1 for row in rows if row["Exposure"] == "Buried")
        variant_row = next((row for row in rows if int(row["Position"]) == variant_position), None)
        points = [
            {
                "position": int(row["Position"]),
                "rsa_height_pct": round(float(row["RSA"]) * 100.0, 2),
                "asa_height_pct": round((float(row["ASA"]) / max_asa) * 100.0, 2),
                "is_variant_site": int(row["Position"]) == variant_position,
            }
            for row in rows
        ]

        return {
            "points": points,
            "mean_rsa": round(mean_rsa, 3),
            "mean_asa": round(mean_asa, 2),
            "exposed_count": exposed_count,
            "buried_count": buried_count,
            "variant_row": variant_row,
            "table_rows": rows[:40],
        }

    @staticmethod
    def _secondary_structure_panel(
        predictor_tables: list[dict[str, Any]], variant_position: int, sequence_length: int
    ) -> dict[str, Any] | None:
        residue_table = next(
            (table for table in predictor_tables if table.get("title") == "PSIPRED residue assignments"),
            None,
        )
        if residue_table is None:
            return None

        raw_rows = residue_table.get("rows", [])
        if not raw_rows:
            return None

        rows = []
        for raw_row in raw_rows:
            state_code = str(raw_row.get("State code", "C")).strip().upper() or "C"
            meta = STRUCTURE_STATE_META.get(state_code, STRUCTURE_STATE_META["C"])
            rows.append(
                {
                    "position": int(raw_row["Position"]),
                    "aa": str(raw_row["AA"]),
                    "state": meta["label"],
                    "state_code": state_code,
                    "class_name": meta["class_name"],
                    "accent": meta["accent"],
                    "confidence": float(raw_row["Confidence"]),
                    "confidence_pct": round(float(raw_row["Confidence"]) * 100.0, 1),
                    "coil_probability": float(raw_row["Coil probability"]),
                    "helix_probability": float(raw_row["Helix probability"]),
                    "strand_probability": float(raw_row["Strand probability"]),
                    "is_variant_site": int(raw_row["Position"]) == variant_position,
                }
            )

        variant_row = next((row for row in rows if row["position"] == variant_position), None)
        if variant_row is None:
            return None

        counts = Counter(row["state_code"] for row in rows)
        total = len(rows)
        mean_confidence = sum(row["confidence"] for row in rows) / total
        state_cards = []
        for state_code in ("H", "E", "C"):
            meta = STRUCTURE_STATE_META[state_code]
            count = counts.get(state_code, 0)
            state_cards.append(
                {
                    "label": meta["label"],
                    "code": state_code,
                    "class_name": meta["class_name"],
                    "count": count,
                    "percent": round((count / total) * 100.0, 1) if total else 0.0,
                }
            )

        segments = SequenceAnalysisService._structure_segments(rows, sequence_length)
        window_radius = 9
        start_index = max(0, variant_position - 1 - window_radius)
        end_index = min(total, variant_position + window_radius)
        neighborhood_rows = rows[start_index:end_index]
        probability_cards = [
            {
                "label": "Coil",
                "class_name": "coil",
                "value": round(variant_row["coil_probability"] * 100.0, 1),
            },
            {
                "label": "Helix",
                "class_name": "helix",
                "value": round(variant_row["helix_probability"] * 100.0, 1),
            },
            {
                "label": "Strand",
                "class_name": "strand",
                "value": round(variant_row["strand_probability"] * 100.0, 1),
            },
        ]
        annotation_blocks = SequenceAnalysisService._structure_annotation_blocks(rows)
        diagram_rows = SequenceAnalysisService._structure_diagram_rows(rows, max_per_line=100)

        return {
            "state_cards": state_cards,
            "variant_row": variant_row,
            "mean_confidence": round(mean_confidence, 3),
            "variant_pct": round((variant_position / sequence_length) * 100.0, 3),
            "segments": segments,
            "neighborhood_rows": neighborhood_rows,
            "probability_cards": probability_cards,
            "annotation_blocks": annotation_blocks,
            "diagram_rows": diagram_rows,
        }

    @staticmethod
    def _structure_segments(rows: list[dict[str, Any]], sequence_length: int) -> list[dict[str, Any]]:
        segments = []
        start_index = 0
        current_state = rows[0]["state_code"]

        for index, row in enumerate(rows[1:], start=1):
            if row["state_code"] == current_state:
                continue
            segments.append(
                SequenceAnalysisService._structure_segment_payload(
                    rows, start_index, index - 1, current_state, sequence_length
                )
            )
            start_index = index
            current_state = row["state_code"]

        segments.append(
            SequenceAnalysisService._structure_segment_payload(
                rows, start_index, len(rows) - 1, current_state, sequence_length
            )
        )
        return segments

    @staticmethod
    def _structure_segment_payload(
        rows: list[dict[str, Any]], start_index: int, end_index: int, state_code: str, sequence_length: int
    ) -> dict[str, Any]:
        meta = STRUCTURE_STATE_META.get(state_code, STRUCTURE_STATE_META["C"])
        start = rows[start_index]["position"]
        end = rows[end_index]["position"]
        return {
            "label": meta["label"],
            "state_code": state_code,
            "class_name": meta["class_name"],
            "start": start,
            "end": end,
            "length": end - start + 1,
            "left_pct": round(((start - 1) / sequence_length) * 100.0, 3),
            "width_pct": max(round(((end - start + 1) / sequence_length) * 100.0, 3), 0.8),
        }

    @staticmethod
    def _structure_annotation_blocks(rows: list[dict[str, Any]], block_size: int = 30) -> list[dict[str, Any]]:
        blocks = []
        for start_index in range(0, len(rows), block_size):
            chunk = rows[start_index : start_index + block_size]
            blocks.append(
                {
                    "start": chunk[0]["position"],
                    "end": chunk[-1]["position"],
                    "residues": chunk,
                }
            )
        return blocks

    @staticmethod
    def _structure_diagram_rows(
        rows: list[dict[str, Any]], max_per_line: int = 100
    ) -> list[dict[str, Any]]:
        """Split residues into lines and group consecutive residues by SS state.

        Each line contains up to *max_per_line* residues.  Within each line
        consecutive residues sharing the same secondary-structure state are
        grouped into a ``segment`` dict so the template can render PDBsum-style
        shapes (helix spring, strand arrow, coil line).
        """
        lines: list[dict[str, Any]] = []
        for line_offset in range(0, len(rows), max_per_line):
            chunk = rows[line_offset: line_offset + max_per_line]
            segments: list[dict[str, Any]] = []
            current_state = chunk[0]["state_code"]
            current_residues: list[dict[str, Any]] = []
            for r in chunk:
                if r["state_code"] == current_state:
                    current_residues.append(r)
                else:
                    meta = STRUCTURE_STATE_META.get(current_state, STRUCTURE_STATE_META["C"])
                    segments.append(
                        {
                            "state_code": current_state,
                            "class_name": meta["class_name"],
                            "label": meta["label"],
                            "residues": current_residues,
                            "length": len(current_residues),
                        }
                    )
                    current_state = r["state_code"]
                    current_residues = [r]
            meta = STRUCTURE_STATE_META.get(current_state, STRUCTURE_STATE_META["C"])
            segments.append(
                {
                    "state_code": current_state,
                    "class_name": meta["class_name"],
                    "label": meta["label"],
                    "residues": current_residues,
                    "length": len(current_residues),
                }
            )
            lines.append(
                {
                    "line_number": line_offset // max_per_line + 1,
                    "start": chunk[0]["position"],
                    "end": chunk[-1]["position"],
                    "segments": segments,
                }
            )
        return lines

    @staticmethod
    def _structural_verdict(
        secondary_structure_panel: dict[str, Any] | None,
        accessibility_panel: dict[str, Any] | None,
        overlapping_intervals: list[dict[str, Any]],
    ) -> dict[str, Any] | None:
        if secondary_structure_panel is None and accessibility_panel is None:
            return None

        parts: list[str] = []
        ss_class: str | None = None
        rsa_class: str | None = None

        if secondary_structure_panel is not None:
            vr = secondary_structure_panel.get("variant_row")
            if vr:
                ss_class = vr.get("class_name")
                conf = vr.get("confidence_pct", 0.0)
                state = vr.get("state", "")
                parts.append(f"{state} (PSIPRED confidence {conf}%)")

        if accessibility_panel is not None:
            vr = accessibility_panel.get("variant_row")
            if vr:
                exposure = str(vr.get("Exposure", ""))
                rsa_class = exposure.lower()
                rsa = float(vr.get("RSA", 0.0))
                parts.append(f"{exposure} (RSA {rsa:.3f})")

        domain_context = ""
        if overlapping_intervals:
            domain_context = overlapping_intervals[0]["label"]

        if not parts:
            return None

        verdict_sentence = "Variant site: " + " · ".join(parts) + "."
        if domain_context:
            verdict_sentence += f" Located in {domain_context}."

        return {
            "parts": parts,
            "domain_context": domain_context,
            "verdict_sentence": verdict_sentence,
            "ss_class": ss_class or "coil",
            "rsa_class": rsa_class or "unknown",
        }

    @staticmethod
    def _amphipathic_helix_panel(
        secondary_structure_panel: dict[str, Any] | None,
        sequence: str,
        variant_position: int,
    ) -> dict[str, Any] | None:
        if secondary_structure_panel is None:
            return None

        segments = secondary_structure_panel.get("segments", [])
        variant_helix = next(
            (
                s
                for s in segments
                if s.get("state_code") == "H"
                and s["start"] <= variant_position <= s["end"]
                and s.get("length", 0) >= 7
            ),
            None,
        )
        if variant_helix is None:
            return None

        seg_start = variant_helix["start"]
        seg_end = variant_helix["end"]
        seg_seq = sequence[seg_start - 1: seg_end]
        mu_h = _hydrophobic_moment(seg_seq)

        if mu_h < 0.35:
            return None

        if mu_h >= 0.60:
            interpretation = "highly amphipathic"
            detail = "consistent with a membrane-inserting or strong protein-interaction helix"
        elif mu_h >= 0.45:
            interpretation = "moderately amphipathic"
            detail = "consistent with a surface-exposed interaction helix"
        else:
            interpretation = "weakly amphipathic"
            detail = "with mild segregation of hydrophobic and hydrophilic faces"

        return {
            "segment_start": seg_start,
            "segment_end": seg_end,
            "segment_length": seg_end - seg_start + 1,
            "hydrophobic_moment": round(mu_h, 3),
            "interpretation": interpretation,
            "summary": (
                f"The helical segment containing the variant (residues {seg_start}–{seg_end}, "
                f"{seg_end - seg_start + 1} aa) is {interpretation} (μH = {mu_h:.2f}), "
                f"{detail}."
            ),
        }

    @staticmethod
    def _supplementary_residue_tracks(residue_tracks: list[dict[str, Any]]) -> list[dict[str, Any]]:
        hidden_keys = {"psipred_confidence", "psipred_helix_probability", "psipred_strand_probability", "dmvfl_rsa", "dmvfl_asa_scaled"}
        return [track for track in residue_tracks if track.get("key") not in hidden_keys]

    @staticmethod
    def _supplementary_predictor_tables(predictor_tables: list[dict[str, Any]]) -> list[dict[str, Any]]:
        hidden_titles = {
            "PSIPRED residue counts",
            "PSIPRED residue assignments",
            "DMVFL-RSA accessibility summary",
            "DMVFL-RSA residue accessibility",
        }
        return [table for table in predictor_tables if table.get("title") not in hidden_titles]


def _charge_class(residue: str) -> str:
    if residue in POSITIVE:
        return "positive"
    if residue in NEGATIVE:
        return "negative"
    return "neutral"


def _polarity_class(residue: str) -> str:
    if residue in POLAR:
        return "polar"
    return "non-polar"


def _hydrophobic_moment(sequence: str, angle_deg: float = 100.0) -> float:
    """Compute the hydrophobic moment (Eisenberg 1982) for an alpha-helix segment.

    angle_deg is the rotation per residue (100° for alpha-helix).
    Returns the normalised hydrophobic moment μH.
    """
    if not sequence:
        return 0.0
    angle_rad = radians(angle_deg)
    sin_sum = sum(HYDROPATHY.get(aa, 0.0) * sin(i * angle_rad) for i, aa in enumerate(sequence))
    cos_sum = sum(HYDROPATHY.get(aa, 0.0) * cos(i * angle_rad) for i, aa in enumerate(sequence))
    return sqrt(sin_sum ** 2 + cos_sum ** 2) / len(sequence)


def _net_charge(sequence: str, ph: float) -> float:
    counts = Counter(sequence)
    positive = (
        1 / (1 + 10 ** (ph - PKA["Nterm"]))
        + counts["K"] / (1 + 10 ** (ph - PKA["K"]))
        + counts["R"] / (1 + 10 ** (ph - PKA["R"]))
        + counts["H"] / (1 + 10 ** (ph - PKA["H"]))
    )
    negative = (
        1 / (1 + 10 ** (PKA["Cterm"] - ph))
        + counts["D"] / (1 + 10 ** (PKA["D"] - ph))
        + counts["E"] / (1 + 10 ** (PKA["E"] - ph))
        + counts["C"] / (1 + 10 ** (PKA["C"] - ph))
        + counts["Y"] / (1 + 10 ** (PKA["Y"] - ph))
    )
    return positive - negative
