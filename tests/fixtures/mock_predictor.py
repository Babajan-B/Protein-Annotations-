from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--variant-position", required=True, type=int)
    args = parser.parse_args()

    fasta = Path(args.input_fasta).read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in fasta if line and not line.startswith(">"))

    payload = {
        "summary": "Mock predictor completed.",
        "interval_features": [
            {
                "label": "Mock domain",
                "category": "domains",
                "start": max(1, args.variant_position - 2),
                "end": min(len(sequence), args.variant_position + 2),
                "variant_hit": True,
            }
        ],
        "site_features": [
            {
                "label": "Mock site",
                "category": "ptm",
                "position": args.variant_position,
                "variant_hit": True,
            }
        ],
        "residue_tracks": [
            {
                "key": "mock_score",
                "label": "Mock score",
                "values": [
                    {"position": index + 1, "value": 0.1 if index + 1 != args.variant_position else 0.9}
                    for index in range(len(sequence))
                ],
            }
        ],
        "evidence_rows": [
            {
                "label": "Mock evidence",
                "position": args.variant_position,
                "summary": "Demonstration evidence row",
            }
        ],
    }
    Path(args.output_json).write_text(json.dumps(payload), encoding="utf-8")


if __name__ == "__main__":
    main()
