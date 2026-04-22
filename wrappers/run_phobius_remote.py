#!/usr/bin/env python3
from __future__ import annotations

import re

from predictor_utils import EbiJobClient, build_parser, interval_feature, read_single_fasta, site_feature, table, write_payload


FT_PATTERN = re.compile(r"^FT\s+(\S+)\s+(\d+)\s+(\d+)\s*(.*?)\s*$")


def main() -> None:
    parser = build_parser("Run Phobius through the official EBI Job Dispatcher service.")
    args = parser.parse_args()

    _, sequence = read_single_fasta(args.input_fasta)
    client = EbiJobClient(
        "phobius",
        request_timeout=args.request_timeout,
        max_wait_seconds=args.max_wait_seconds,
        poll_interval_seconds=args.poll_interval_seconds,
    )
    job_id = client.submit({"sequence": f">query\n{sequence}\n", "format": "long"})
    client.wait(job_id)
    raw_output = client.result(job_id, "out")
    intervals, sites = _parse_output(raw_output, args.variant_position)

    payload = {
        "summary": f"Phobius completed with {len(intervals)} topology intervals.",
        "interval_features": intervals,
        "site_features": sites,
        "evidence_rows": [
            {
                "source": "Phobius",
                "kind": "topology",
                "label": "Topology summary",
                "position": args.variant_position,
                "summary": _variant_summary(intervals, sites, args.variant_position),
                "score": str(len(intervals)),
            }
        ],
        "tables": [
            table(
                "Phobius topology features",
                ["Type", "Start", "End", "Description"],
                [
                    {
                        "Type": feature["label"],
                        "Start": feature["start"],
                        "End": feature["end"],
                        "Description": feature["description"],
                    }
                    for feature in intervals
                ],
            )
        ],
        "notices": client.default_notices(),
    }
    write_payload(args.output_json, payload)


def _parse_output(raw_output: str, variant_position: int) -> tuple[list[dict], list[dict]]:
    intervals: list[dict] = []
    sites: list[dict] = []

    for raw_line in raw_output.splitlines():
        match = FT_PATTERN.match(raw_line)
        if not match:
            continue
        feature_type, start_str, end_str, description = match.groups()
        start = int(start_str)
        end = int(end_str)
        cleaned_description = description.strip().strip(".")

        if feature_type == "SIGNAL":
            intervals.append(
                interval_feature(
                    source="Phobius",
                    label="Signal peptide",
                    category="topology",
                    start=start,
                    end=end,
                    variant_position=variant_position,
                    type_name="SIGNAL_PEPTIDE",
                    description="Predicted N-terminal signal peptide.",
                )
            )
            sites.append(
                site_feature(
                    source="Phobius",
                    label="Signal peptide cleavage site",
                    category="topology",
                    position=end,
                    variant_position=variant_position,
                    type_name="SIGNAL_CLEAVAGE",
                    description=f"Predicted signal peptide ends at residue {end}.",
                )
            )
            continue

        label = _label_for_feature(feature_type, cleaned_description)
        intervals.append(
            interval_feature(
                source="Phobius",
                label=label,
                category="topology",
                start=start,
                end=end,
                variant_position=variant_position,
                type_name=f"PHOBIUS_{feature_type}",
                description=cleaned_description or label,
            )
        )

    return intervals, sites


def _label_for_feature(feature_type: str, description: str) -> str:
    if feature_type == "TRANSMEM":
        return "Transmembrane helix"
    if feature_type == "DOMAIN" and description:
        return description.title()
    return feature_type.title()


def _variant_summary(intervals: list[dict], sites: list[dict], variant_position: int) -> str:
    hits = [feature["label"] for feature in intervals if feature["start"] <= variant_position <= feature["end"]]
    if any(site["position"] == variant_position for site in sites):
        hits.append("cleavage site")
    if hits:
        return f"Residue {variant_position} overlaps: {', '.join(hits)}."
    return f"Residue {variant_position} does not overlap a predicted Phobius feature."


if __name__ == "__main__":
    main()
