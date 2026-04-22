from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
from io import StringIO
from typing import Any


def serialize_report(report: dict[str, Any]) -> dict[str, Any]:
    return _prune_empty_tool_sections(_convert(report))


def report_to_csv(report: dict[str, Any]) -> str:
    output = StringIO()
    writer = csv.writer(output)

    writer.writerow(["section", "field", "value"])
    writer.writerow(["protein", "accession", report["protein"].accession])
    writer.writerow(["protein", "entry_name", report["protein"].entry_name])
    writer.writerow(["protein", "gene_symbol", report["protein"].gene_symbol])
    writer.writerow(["protein", "protein_name", report["protein"].protein_name])
    writer.writerow(["protein", "organism_name", report["protein"].organism_name])
    writer.writerow(["protein", "length", report["protein"].length])
    writer.writerow(["variant", "submitted", report["variant"].original_text])
    writer.writerow(["variant", "normalized", report["variant"].short_notation])
    writer.writerow(["variant", "hgvs", report["variant"].hgvs_protein])

    for metric in report["metric_deltas"]:
        writer.writerow(["metric", metric.label, metric.delta])

    for feature in report["interval_features"]:
        writer.writerow(
            [
                "interval_feature",
                feature["label"],
                f'{feature["source"]}|{feature["type"]}|{feature["start"]}-{feature["end"]}',
            ]
        )

    for feature in report["site_features"]:
        writer.writerow(
            [
                "site_feature",
                feature["label"],
                f'{feature["source"]}|{feature["type"]}|{feature["position"]}',
            ]
        )

    for evidence in report["evidence_rows"]:
        writer.writerow(
            [
                "evidence",
                evidence["label"],
                f'{evidence["source"]}|{evidence["position"]}|{evidence["summary"]}',
            ]
        )

    return output.getvalue()


def structure_annotations_to_csv(report: dict[str, Any]) -> str:
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(
        [
            "position",
            "amino_acid",
            "state",
            "state_code",
            "confidence",
            "coil_probability",
            "helix_probability",
            "strand_probability",
            "is_variant_site",
        ]
    )

    panel = report.get("secondary_structure_panel") or {}
    for block in panel.get("annotation_blocks", []):
        for residue in block.get("residues", []):
            writer.writerow(
                [
                    residue["position"],
                    residue["aa"],
                    residue["state"],
                    residue["state_code"],
                    residue["confidence"],
                    residue["coil_probability"],
                    residue["helix_probability"],
                    residue["strand_probability"],
                    "yes" if residue["is_variant_site"] else "no",
                ]
            )

    return output.getvalue()


def _convert(value: Any) -> Any:
    if is_dataclass(value):
        return _convert(asdict(value))
    if isinstance(value, dict):
        return {key: _convert(inner) for key, inner in value.items()}
    if isinstance(value, list):
        return [_convert(item) for item in value]
    return value


def _prune_empty_tool_sections(report: dict[str, Any]) -> dict[str, Any]:
    removable_keys = {
        "accessibility_panel",
        "architecture_rows",
        "chromosome_location_panel",
        "conservation_panel",
        "evidence_plot",
        "evidence_rows",
        "feature_counts",
        "insight_notices",
        "interval_features",
        "notices",
        "phylogenetic_panel",
        "predictor_tables",
        "residue_tracks",
        "resource_cards",
        "secondary_structure_panel",
        "site_features",
        "source_status",
        "structure_architecture_rows",
        "structure_interval_features",
        "structure_site_features",
        "supplementary_predictor_tables",
        "supplementary_residue_tracks",
    }
    pruned = {}
    for key, value in report.items():
        if key in removable_keys and (value is None or value == [] or value == {}):
            continue
        pruned[key] = value
    return pruned
