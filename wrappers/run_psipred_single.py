#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from pathlib import Path
import subprocess
import tempfile

from predictor_utils import build_parser, interval_feature, read_single_fasta, site_feature, table, write_payload


STRUCTURE_LABELS = {"H": "Alpha helix", "E": "Beta strand", "C": "Coil"}


def main() -> None:
    parser = build_parser("Run PSIPRED single-sequence secondary-structure prediction.")
    parser.add_argument("--psipred-root", default="")
    args = parser.parse_args()

    _, sequence = read_single_fasta(args.input_fasta)
    repo_root = Path(__file__).resolve().parents[1]
    psipred_root = Path(args.psipred_root) if args.psipred_root else repo_root / "envs" / "bio-tools"
    bin_dir = psipred_root / "bin"
    data_dir = psipred_root / "share" / "psipred" / "data"
    cache_dir = repo_root / ".analysis-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="psipred-", dir=cache_dir) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        fasta_path = temp_dir / "input.fasta"
        mtx_path = temp_dir / "input.mtx"
        ss_path = temp_dir / "out.ss"
        ss2_path = temp_dir / "out.ss2"
        horiz_path = temp_dir / "out.horiz"
        fasta_path.write_text(f">query\n{sequence}\n", encoding="utf-8")

        with mtx_path.open("w", encoding="utf-8") as handle:
            subprocess.run([str(bin_dir / "seq2mtx"), str(fasta_path)], check=True, stdout=handle)
        with ss_path.open("w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    str(bin_dir / "psipred"),
                    str(mtx_path),
                    str(data_dir / "weights.dat"),
                    str(data_dir / "weights.dat2"),
                    str(data_dir / "weights.dat3"),
                ],
                check=True,
                stdout=handle,
            )
        with horiz_path.open("w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    str(bin_dir / "psipass2"),
                    str(data_dir / "weights_p2.dat"),
                    "1",
                    "1.0",
                    "1.0",
                    str(ss2_path),
                    str(ss_path),
                ],
                check=True,
                stdout=handle,
            )

        rows = _parse_ss2(ss2_path)

    interval_features = _build_interval_features(rows, args.variant_position)
    variant_row = rows[args.variant_position - 1]
    counts = Counter(row["prediction"] for row in rows)

    payload = {
        "summary": (
            f"PSIPRED single-sequence run completed with {counts['H']} helix residues and "
            f"{counts['E']} strand residues."
        ),
        "interval_features": interval_features,
        "site_features": [
            site_feature(
                source="PSIPRED",
                label=f"Predicted {STRUCTURE_LABELS[variant_row['prediction']]}",
                category="structure",
                position=args.variant_position,
                variant_position=args.variant_position,
                type_name=f"PSIPRED_{variant_row['prediction']}",
                description=(
                    f"Predicted state at residue {args.variant_position} is "
                    f"{STRUCTURE_LABELS[variant_row['prediction']].lower()} "
                    f"with confidence {variant_row['confidence']:.2f}."
                ),
                evidence=f"PSIPRED probabilities C={variant_row['coil']:.3f}, H={variant_row['helix']:.3f}, E={variant_row['strand']:.3f}",
            )
        ],
        "residue_tracks": [
            {
                "key": "psipred_confidence",
                "label": "PSIPRED confidence",
                "values": [{"position": row["position"], "value": row["confidence"]} for row in rows],
            },
            {
                "key": "psipred_helix_probability",
                "label": "Helix probability",
                "values": [{"position": row["position"], "value": row["helix"]} for row in rows],
            },
            {
                "key": "psipred_strand_probability",
                "label": "Strand probability",
                "values": [{"position": row["position"], "value": row["strand"]} for row in rows],
            },
        ],
        "evidence_rows": [
            {
                "source": "PSIPRED",
                "kind": "secondary_structure",
                "label": "Variant-site structural context",
                "position": args.variant_position,
                "summary": (
                    f"Residue {args.variant_position} is predicted as "
                    f"{STRUCTURE_LABELS[variant_row['prediction']].lower()} "
                    f"with confidence {variant_row['confidence']:.2f}."
                ),
                "score": f"{variant_row['confidence']:.3f}",
            }
        ],
        "tables": [
            table(
                "PSIPRED residue counts",
                ["State", "Residues"],
                [
                    {"State": "Alpha helix", "Residues": counts["H"]},
                    {"State": "Beta strand", "Residues": counts["E"]},
                    {"State": "Coil", "Residues": counts["C"]},
                ],
            ),
            table(
                "PSIPRED residue assignments",
                [
                    "Position",
                    "AA",
                    "State",
                    "State code",
                    "Confidence",
                    "Coil probability",
                    "Helix probability",
                    "Strand probability",
                ],
                [
                    {
                        "Position": row["position"],
                        "AA": row["residue"],
                        "State": STRUCTURE_LABELS[row["prediction"]],
                        "State code": row["prediction"],
                        "Confidence": round(float(row["confidence"]), 4),
                        "Coil probability": round(float(row["coil"]), 4),
                        "Helix probability": round(float(row["helix"]), 4),
                        "Strand probability": round(float(row["strand"]), 4),
                    }
                    for row in rows
                ],
            ),
        ],
    }
    write_payload(args.output_json, payload)


def _parse_ss2(path: Path) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        position_str, residue, prediction, coil, helix, strand = line.split()
        probabilities = [float(coil), float(helix), float(strand)]
        rows.append(
            {
                "position": int(position_str),
                "residue": residue,
                "prediction": prediction,
                "coil": probabilities[0],
                "helix": probabilities[1],
                "strand": probabilities[2],
                "confidence": max(probabilities),
            }
        )
    return rows


def _build_interval_features(rows: list[dict[str, float | int | str]], variant_position: int) -> list[dict[str, object]]:
    features: list[dict[str, object]] = []
    segment_start = 1
    current_state = rows[0]["prediction"]

    for index, row in enumerate(rows[1:], start=2):
        if row["prediction"] == current_state:
            continue
        features.extend(_segment_feature(rows, segment_start, index - 1, current_state, variant_position))
        segment_start = index
        current_state = row["prediction"]

    features.extend(_segment_feature(rows, segment_start, len(rows), current_state, variant_position))
    return features


def _segment_feature(
    rows: list[dict[str, float | int | str]],
    start: int,
    end: int,
    state: str,
    variant_position: int,
) -> list[dict[str, object]]:
    if state == "C" or (end - start + 1) < 2:
        return []
    confidence = sum(float(rows[index - 1]["confidence"]) for index in range(start, end + 1)) / (end - start + 1)
    label = STRUCTURE_LABELS[str(state)]
    return [
        interval_feature(
            source="PSIPRED",
            label=label,
            category="structure",
            start=start,
            end=end,
            variant_position=variant_position,
            type_name=f"PSIPRED_{state}",
            description=f"{label} segment predicted by PSIPRED with mean confidence {confidence:.2f}.",
            evidence=f"Mean confidence {confidence:.3f}",
        )
    ]


if __name__ == "__main__":
    main()
