#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from io import StringIO
import re

from Bio import AlignIO

from predictor_utils import AMINO_ACIDS, EbiJobClient, build_parser, read_single_fasta, site_feature, table, write_payload


HIT_LINE_PATTERN = re.compile(
    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(.*)$"
)


def main() -> None:
    parser = build_parser("Run HMMER phmmer against SwissProt and derive per-residue conservation.")
    parser.add_argument("--database", default="swissprot")
    parser.add_argument("--nhits", type=int, default=64)
    args = parser.parse_args()

    _, sequence = read_single_fasta(args.input_fasta)
    client = EbiJobClient(
        "hmmer3_phmmer",
        request_timeout=args.request_timeout,
        max_wait_seconds=args.max_wait_seconds,
        poll_interval_seconds=args.poll_interval_seconds,
    )
    job_id = client.submit(
        {
            "sequence": f">query\n{sequence}\n",
            "database": args.database,
            "nhits": str(args.nhits),
        }
    )
    client.wait(job_id)
    stockholm = client.result(job_id, "sto")
    out_text = client.result(job_id, "out")
    track, informative_depth = _build_conservation_track(stockholm)
    variant_point = track[args.variant_position - 1]
    top_hits = _parse_hits(out_text)

    if variant_point["value"] >= 0.8:
        label = "Highly conserved wild-type residue"
    elif variant_point["value"] >= 0.5:
        label = "Moderately conserved wild-type residue"
    else:
        label = "Variable residue across homologs"

    payload = {
        "summary": (
            f"phmmer recovered {len(top_hits)} homologs and computed conservation over "
            f"{informative_depth} aligned sequences."
        ),
        "site_features": [
            site_feature(
                source="HMMER phmmer",
                label=label,
                category="conservation",
                position=args.variant_position,
                variant_position=args.variant_position,
                type_name="CONSERVATION_SITE",
                description=(
                    f"Wild-type residue frequency at site {args.variant_position}: "
                    f"{variant_point['value']:.2%} across {variant_point['depth']} informative homologs."
                ),
                evidence=f"Consensus residue {variant_point['consensus']}, depth {variant_point['depth']}",
            )
        ],
        "residue_tracks": [
            {
                "key": "phmmer_match_fraction",
                "label": "Homolog match fraction",
                "values": [{"position": point["position"], "value": point["value"]} for point in track],
            }
        ],
        "evidence_rows": [
            {
                "source": "HMMER phmmer",
                "kind": "conservation",
                "label": "Variant-site conservation",
                "position": args.variant_position,
                "summary": (
                    f"Wild-type residue {sequence[args.variant_position - 1]} is observed in "
                    f"{variant_point['value']:.2%} of {variant_point['depth']} informative homologs."
                ),
                "score": f"{variant_point['value']:.3f}",
            }
        ],
        "tables": [
            table(
                "Top phmmer homologs",
                ["Accession", "E-value", "Score", "Description"],
                top_hits[:15],
            ),
            table(
                "Conservation summary",
                ["Metric", "Value"],
                [
                    {"Metric": "Alignment depth", "Value": informative_depth},
                    {"Metric": "Variant-site match fraction", "Value": f"{variant_point['value']:.2%}"},
                    {"Metric": "Variant-site consensus", "Value": variant_point["consensus"]},
                ],
            ),
        ],
        "notices": client.default_notices(),
    }
    write_payload(args.output_json, payload)


def _build_conservation_track(stockholm_text: str) -> tuple[list[dict[str, float | int | str]], int]:
    alignment = AlignIO.read(StringIO(stockholm_text), "stockholm")
    query_record = alignment[0]
    homolog_records = alignment[1:] if len(alignment) > 1 else alignment
    track: list[dict[str, float | int | str]] = []
    query_position = 0
    max_depth = 0

    for column_index in range(alignment.get_alignment_length()):
        query_residue = query_record.seq[column_index].upper()
        if query_residue not in AMINO_ACIDS:
            continue
        query_position += 1
        column_residues = [
            record.seq[column_index].upper()
            for record in homolog_records
            if record.seq[column_index].upper() in AMINO_ACIDS
        ]
        if not column_residues:
            match_fraction = 0.0
            consensus = query_residue
        else:
            counts = Counter(column_residues)
            consensus, _ = counts.most_common(1)[0]
            match_fraction = counts.get(query_residue, 0) / len(column_residues)
            max_depth = max(max_depth, len(column_residues))
        track.append(
            {
                "position": query_position,
                "value": round(match_fraction, 4),
                "consensus": consensus,
                "depth": len(column_residues),
            }
        )

    return track, max_depth


def _parse_hits(out_text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    table_started = False
    for raw_line in out_text.splitlines():
        if "Scores for complete sequences" in raw_line:
            table_started = True
            continue
        if not table_started or not raw_line.strip():
            continue
        if raw_line.strip().startswith("Domain annotation"):
            break
        match = HIT_LINE_PATTERN.match(raw_line)
        if not match:
            continue
        full_evalue, score, _bias, _best_evalue, _best_score, _best_bias, _exp, _n, accession, description = match.groups()
        rows.append(
            {
                "Accession": accession,
                "E-value": full_evalue,
                "Score": score,
                "Description": description,
            }
        )
    return rows


if __name__ == "__main__":
    main()
