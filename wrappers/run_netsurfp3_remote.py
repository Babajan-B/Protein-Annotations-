#!/usr/bin/env python3
"""
NetSurfP-3.0 wrapper — secondary structure, RSA, and disorder from local install.

Install NetSurfP-3.0 (https://services.healthtech.dtu.dk/services/NetSurfP-3.0/)
then set NETSURFP3_COMMAND_TEMPLATE to point at this wrapper:

  NETSURFP3_COMMAND_TEMPLATE=python /path/to/run_netsurfp3_remote.py \\
      --input-fasta {input_fasta} --output-json {output_json} \\
      --variant-position {variant_position}

The netsurfp3 binary is searched in PATH unless overridden with --netsurfp3-bin.
"""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from predictor_utils import (
    build_parser,
    interval_feature,
    read_single_fasta,
    table,
    write_payload,
)


SS_LABELS = {"H": "Alpha helix", "E": "Beta strand", "C": "Coil"}


def main() -> None:
    parser = build_parser(
        "Run NetSurfP-3.0 for secondary structure, RSA, and disorder prediction."
    )
    parser.add_argument(
        "--netsurfp3-bin",
        default="netsurfp3",
        help="Path to the netsurfp3 executable (default: netsurfp3 in PATH).",
    )
    args = parser.parse_args()

    _, _sequence = read_single_fasta(args.input_fasta)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_csv = Path(tmpdir) / "netsurfp3_output.csv"
        cmd = [
            args.netsurfp3_bin,
            "--csv",
            str(out_csv),
            args.input_fasta,
        ]
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "netsurfp3 binary not found. Install NetSurfP-3.0 and set "
                "--netsurfp3-bin or add it to PATH."
            )
        if proc.returncode != 0:
            raise RuntimeError(
                f"netsurfp3 exited with code {proc.returncode}: "
                f"{(proc.stderr or proc.stdout)[:400].strip()}"
            )

        if not out_csv.exists():
            raise RuntimeError("netsurfp3 did not produce the expected CSV output.")

        raw_rows = list(csv.DictReader(out_csv.read_text(encoding="utf-8")))

    if not raw_rows:
        raise RuntimeError("netsurfp3 output CSV is empty.")

    rows = [_normalise_row(r) for r in raw_rows]
    interval_feats = _build_ss_intervals(rows, args.variant_position)
    ss_table_rows, rsa_track, disorder_track = _build_tracks(rows)
    variant_summary = _variant_summary(rows, args.variant_position)

    payload = {
        "summary": f"NetSurfP-3.0 completed: {len(rows)} residues annotated. {variant_summary}",
        "interval_features": interval_feats,
        "site_features": [],
        "residue_tracks": [
            {"key": "netsurfp3_rsa", "label": "NetSurfP-3 RSA", "values": rsa_track},
            {
                "key": "netsurfp3_disorder",
                "label": "NetSurfP-3 Disorder",
                "values": disorder_track,
            },
        ],
        "evidence_rows": [
            {
                "source": "NetSurfP-3.0",
                "kind": "structure",
                "label": "SS / RSA / disorder summary",
                "position": args.variant_position,
                "summary": variant_summary,
                "score": "",
            }
        ],
        "tables": [
            table(
                "NetSurfP-3 residue assignments",
                ["Position", "AA", "SS", "RSA", "Disorder"],
                ss_table_rows,
            )
        ],
        "notices": [
            "NetSurfP-3.0 provides per-residue secondary structure (H/E/C), "
            "relative solvent accessibility, backbone torsion angles, and disorder probability."
        ],
    }
    write_payload(args.output_json, payload)


# ------------------------------------------------------------------ #
#  Row normalisation (handles variations in column naming)             #
# ------------------------------------------------------------------ #


def _normalise_row(row: dict[str, str]) -> dict[str, str | int | float]:
    # Support both DTU CSV column names and common variants
    pos = int(row.get("n", row.get("position", row.get("pos", 0))))
    aa = row.get("aa", row.get("residue", "X"))
    ss = row.get("q3", row.get("sec_str", row.get("ss", "C"))).upper() or "C"
    rsa = float(row.get("rsa", row.get("asa_norm", row.get("rel_asa", 0.0))))
    disorder = float(row.get("disorder", row.get("disor", 0.0)))
    return {"position": pos, "aa": aa, "ss": ss, "rsa": rsa, "disorder": disorder}


# ------------------------------------------------------------------ #
#  Track and table builders                                            #
# ------------------------------------------------------------------ #


def _build_tracks(
    rows: list[dict],
) -> tuple[list[dict], list[dict], list[dict]]:
    ss_table: list[dict] = []
    rsa_vals: list[dict] = []
    disorder_vals: list[dict] = []
    for r in rows:
        ss_table.append(
            {
                "Position": r["position"],
                "AA": r["aa"],
                "SS": SS_LABELS.get(str(r["ss"]), str(r["ss"])),
                "RSA": round(float(r["rsa"]), 3),
                "Disorder": round(float(r["disorder"]), 3),
            }
        )
        rsa_vals.append({"position": r["position"], "value": float(r["rsa"])})
        disorder_vals.append(
            {"position": r["position"], "value": float(r["disorder"])}
        )
    return ss_table, rsa_vals, disorder_vals


def _build_ss_intervals(rows: list[dict], variant_position: int) -> list[dict]:
    if not rows:
        return []
    intervals: list[dict] = []
    current_ss = str(rows[0]["ss"]).upper()
    seg_start = int(rows[0]["position"])
    seg_end = seg_start

    for r in rows[1:]:
        ss = str(r["ss"]).upper()
        pos = int(r["position"])
        if ss == current_ss:
            seg_end = pos
        else:
            if current_ss in ("H", "E"):
                intervals.append(
                    interval_feature(
                        source="NetSurfP-3.0",
                        label=SS_LABELS.get(current_ss, current_ss),
                        category="structure",
                        start=seg_start,
                        end=seg_end,
                        variant_position=variant_position,
                        type_name="secondary_structure",
                    )
                )
            current_ss = ss
            seg_start = seg_end = pos

    if current_ss in ("H", "E"):
        intervals.append(
            interval_feature(
                source="NetSurfP-3.0",
                label=SS_LABELS.get(current_ss, current_ss),
                category="structure",
                start=seg_start,
                end=seg_end,
                variant_position=variant_position,
                type_name="secondary_structure",
            )
        )
    return intervals


def _variant_summary(rows: list[dict], variant_position: int) -> str:
    for r in rows:
        if int(r["position"]) == variant_position:
            ss = str(r["ss"]).upper()
            rsa = float(r["rsa"])
            disorder = float(r["disorder"])
            exposure = "Exposed" if rsa > 0.25 else "Buried"
            return (
                f"Position {variant_position}: {SS_LABELS.get(ss, ss)}, "
                f"RSA {rsa:.3f} ({exposure}), disorder {disorder:.3f}."
            )
    return f"Position {variant_position} not found in NetSurfP-3.0 output."


if __name__ == "__main__":
    main()
