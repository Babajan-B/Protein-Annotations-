#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter

from predictor_utils import EbiJobClient, build_parser, interval_feature, read_single_fasta, site_feature, table, write_payload


def main() -> None:
    parser = build_parser("Run InterProScan 5 through the official EBI Job Dispatcher service.")
    args = parser.parse_args()

    _, sequence = read_single_fasta(args.input_fasta)
    client = EbiJobClient(
        "iprscan5",
        request_timeout=args.request_timeout,
        max_wait_seconds=args.max_wait_seconds,
        poll_interval_seconds=args.poll_interval_seconds,
    )
    job_id = client.submit(
        {
            "sequence": f">query\n{sequence}\n",
            "goterms": "false",
            "pathways": "false",
            "stype": "p",
        }
    )
    client.wait(job_id)
    result_types = client.result_types(job_id)
    if "tsv" not in result_types:
        raise RuntimeError(f"iprscan5 did not provide TSV output. Available result types: {result_types}")
    tsv_output = client.result(job_id, "tsv")

    interval_features, overlapping_rows = _parse_tsv(tsv_output, args.variant_position)
    category_counts = Counter(feature["category"] for feature in interval_features)

    payload = {
        "summary": f"InterProScan returned {len(interval_features)} annotated intervals across {len(category_counts)} categories.",
        "interval_features": interval_features,
        "site_features": [
            site_feature(
                source="InterProScan",
                label="Variant overlaps annotated signature"
                if overlapping_rows
                else "No annotated signature overlap",
                category="domains",
                position=args.variant_position,
                variant_position=args.variant_position,
                type_name="INTERPROSCAN_SITE",
                description=_variant_summary(overlapping_rows, args.variant_position),
            )
        ],
        "evidence_rows": [
            {
                "source": "InterProScan",
                "kind": "signature_overlap",
                "label": "InterProScan overlap count",
                "position": args.variant_position,
                "summary": _variant_summary(overlapping_rows, args.variant_position),
                "score": str(len(overlapping_rows)),
            }
        ],
        "tables": [
            table(
                "InterProScan annotations",
                ["Analysis", "Signature", "Label", "Start", "End", "Category"],
                overlapping_rows[:20] if overlapping_rows else _table_rows(interval_features[:20]),
            )
        ],
        "notices": client.default_notices(),
    }
    write_payload(args.output_json, payload)


def _parse_tsv(tsv_output: str, variant_position: int) -> tuple[list[dict], list[dict]]:
    interval_features: list[dict] = []
    overlapping_rows: list[dict] = []
    seen: set[tuple[str, str, int, int, str]] = set()

    for raw_line in tsv_output.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 13:
            continue
        analysis = columns[3]
        signature_accession = columns[4]
        signature_description = columns[5] or signature_accession
        start = int(columns[6])
        end = int(columns[7])
        interpro_accession = columns[11]
        interpro_description = columns[12]
        label = interpro_description or signature_description or signature_accession
        category = _category_for_analysis(analysis)
        key = (analysis, signature_accession, start, end, label)
        if key in seen:
            continue
        seen.add(key)
        feature = interval_feature(
            source=f"InterProScan:{analysis}",
            label=label,
            category=category,
            start=start,
            end=end,
            variant_position=variant_position,
            type_name=signature_accession or analysis,
            description=signature_description,
            evidence=interpro_accession,
        )
        interval_features.append(feature)
        if start <= variant_position <= end:
            overlapping_rows.append(
                {
                    "Analysis": analysis,
                    "Signature": signature_accession,
                    "Label": label,
                    "Start": start,
                    "End": end,
                    "Category": category,
                }
            )

    return interval_features, overlapping_rows


def _category_for_analysis(analysis: str) -> str:
    lowered = analysis.lower()
    if lowered.startswith("signalp") or lowered == "phobius":
        return "topology"
    if lowered.startswith("mobidblite"):
        return "disorder"
    return "domains"


def _variant_summary(overlapping_rows: list[dict], variant_position: int) -> str:
    if not overlapping_rows:
        return f"Residue {variant_position} does not overlap an InterProScan signature hit."
    labels = ", ".join(row["Label"] for row in overlapping_rows[:4])
    return f"Residue {variant_position} overlaps {len(overlapping_rows)} InterProScan hit(s): {labels}."


def _table_rows(features: list[dict]) -> list[dict]:
    return [
        {
            "Analysis": feature["source"].split(":", 1)[-1],
            "Signature": feature["type"],
            "Label": feature["label"],
            "Start": feature["start"],
            "End": feature["end"],
            "Category": feature["category"],
        }
        for feature in features
    ]


if __name__ == "__main__":
    main()
