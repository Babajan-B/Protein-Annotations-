#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from io import StringIO
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import types

from Bio import AlignIO
import numpy as np
import torch

from predictor_utils import AMINO_ACIDS, EbiJobClient, build_parser, interval_feature, read_single_fasta, site_feature, table, write_payload


AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")
AA_INDEX = {residue: index for index, residue in enumerate(AA_ORDER)}
BACKGROUND = {
    "A": 0.078,
    "R": 0.051,
    "N": 0.043,
    "D": 0.052,
    "C": 0.024,
    "Q": 0.034,
    "E": 0.063,
    "G": 0.073,
    "H": 0.026,
    "I": 0.052,
    "L": 0.091,
    "K": 0.058,
    "M": 0.024,
    "F": 0.040,
    "P": 0.051,
    "S": 0.071,
    "T": 0.058,
    "W": 0.014,
    "Y": 0.033,
    "V": 0.065,
}
ASA_MAX = {
    "A": 115.0,
    "C": 135.0,
    "D": 150.0,
    "E": 190.0,
    "F": 210.0,
    "G": 75.0,
    "H": 195.0,
    "I": 175.0,
    "K": 200.0,
    "L": 170.0,
    "M": 185.0,
    "N": 160.0,
    "P": 145.0,
    "Q": 180.0,
    "R": 225.0,
    "S": 115.0,
    "T": 140.0,
    "V": 155.0,
    "W": 255.0,
    "Y": 230.0,
}
RSA_EXPOSED_THRESHOLD = 0.25
RSA_BURIED_THRESHOLD = 0.1
DEFAULT_LOCAL_STEP_TIMEOUT_SECONDS = 90
DEFAULT_MAX_SEQUENCE_LENGTH = 700


def main() -> None:
    parser = build_parser("Run the local DMVFL-RSA model for residue-wise RSA and ASA prediction.")
    parser.add_argument("--dmvfl-root", default="")
    parser.add_argument("--nhits", type=int, default=128)
    parser.add_argument("--local-step-timeout-seconds", type=int, default=DEFAULT_LOCAL_STEP_TIMEOUT_SECONDS)
    parser.add_argument("--max-sequence-length", type=int, default=DEFAULT_MAX_SEQUENCE_LENGTH)
    args = parser.parse_args()

    _, sequence = read_single_fasta(args.input_fasta)
    if len(sequence) > args.max_sequence_length:
        raise RuntimeError(
            "DMVFL-RSA is disabled for sequences longer than "
            f"{args.max_sequence_length} aa in this app build; received {len(sequence)} aa."
        )
    repo_root = Path(__file__).resolve().parents[1]
    dmvfl_root = Path(args.dmvfl_root) if args.dmvfl_root else repo_root / "DMVFL-RSA-main"
    bio_root = repo_root / "envs" / "bio-tools"

    stockholm = _run_phmmer(sequence, args)
    psfm_matrix = _psfm_from_stockholm(stockholm, sequence)
    pssm_matrix = _pssm_from_psfm(psfm_matrix)
    pss_matrix = _run_psipred(sequence, repo_root, bio_root, args.local_step_timeout_seconds)
    jpsfm_matrix = _run_threader(sequence, psfm_matrix, repo_root, dmvfl_root, args.local_step_timeout_seconds)
    rsa_values = _run_model(pssm_matrix, psfm_matrix[:, :20], pss_matrix, jpsfm_matrix, dmvfl_root)
    asa_values = np.array([ASA_MAX.get(residue, 180.0) * rsa for residue, rsa in zip(sequence, rsa_values)], dtype=float)

    payload = _build_payload(sequence, rsa_values, asa_values, args.variant_position, args)
    write_payload(args.output_json, payload)


def _run_phmmer(sequence: str, args) -> str:
    client = EbiJobClient(
        "hmmer3_phmmer",
        request_timeout=args.request_timeout,
        max_wait_seconds=args.max_wait_seconds,
        poll_interval_seconds=args.poll_interval_seconds,
    )
    job_id = client.submit(
        {
            "sequence": f">query\n{sequence}\n",
            "database": "swissprot",
            "nhits": str(args.nhits),
        }
    )
    client.wait(job_id)
    return client.result(job_id, "sto")


def _psfm_from_stockholm(stockholm_text: str, sequence: str) -> np.ndarray:
    alignment = AlignIO.read(StringIO(stockholm_text), "stockholm")
    query_record = alignment[0]
    homologs = alignment[1:] if len(alignment) > 1 else alignment[0:1]
    rows: list[list[float]] = []

    for column_index in range(alignment.get_alignment_length()):
        query_residue = str(query_record.seq[column_index]).upper()
        if query_residue not in AMINO_ACIDS:
            continue

        residue_counts = Counter()
        gap_count = 0
        for record in homologs:
            residue = str(record.seq[column_index]).upper()
            if residue in AMINO_ACIDS:
                residue_counts[residue] += 1
            else:
                gap_count += 1

        if not residue_counts:
            row = [1.0 if residue == query_residue else 0.0 for residue in AA_ORDER]
            row.append(0.0)
            rows.append(row)
            continue

        depth = sum(residue_counts.values())
        row = [residue_counts.get(residue, 0) / depth for residue in AA_ORDER]
        row.append(gap_count / max(len(homologs), 1))
        rows.append(row)

    matrix = np.array(rows, dtype=np.float32)
    if matrix.shape[0] != len(sequence):
        raise RuntimeError(f"Alignment-derived PSFM length {matrix.shape[0]} does not match sequence length {len(sequence)}.")
    return matrix


def _pssm_from_psfm(psfm_matrix: np.ndarray) -> np.ndarray:
    pssm = np.zeros((psfm_matrix.shape[0], 20), dtype=np.float32)
    for index, row in enumerate(psfm_matrix[:, :20]):
        for aa_index, residue in enumerate(AA_ORDER):
            freq = max(float(row[aa_index]), 1e-5)
            background = BACKGROUND[residue]
            score = math.log(freq / background)
            pssm[index, aa_index] = 1.0 / (1.0 + math.exp(-(score * 2.0)))
    return pssm


def _run_psipred(sequence: str, repo_root: Path, bio_root: Path, step_timeout_seconds: int) -> np.ndarray:
    cache_dir = repo_root / ".analysis-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    bin_dir = bio_root / "bin"
    data_dir = bio_root / "share" / "psipred" / "data"
    with tempfile.TemporaryDirectory(prefix="dmvfl-psipred-", dir=cache_dir) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        fasta_path = temp_dir / "input.fasta"
        mtx_path = temp_dir / "input.mtx"
        ss_path = temp_dir / "out.ss"
        ss2_path = temp_dir / "out.ss2"
        horiz_path = temp_dir / "out.horiz"
        fasta_path.write_text(f">query\n{sequence}\n", encoding="utf-8")

        with mtx_path.open("w", encoding="utf-8") as handle:
            subprocess.run(
                [str(bin_dir / "seq2mtx"), str(fasta_path)],
                check=True,
                stdout=handle,
                timeout=step_timeout_seconds,
            )
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
                timeout=step_timeout_seconds,
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
                timeout=step_timeout_seconds,
            )

        rows = []
        for raw_line in ss2_path.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            _position, _residue, _pred, coil, helix, strand = line.split()
            rows.append([float(coil), float(helix), float(strand)])
    return np.array(rows, dtype=np.float32)


def _run_threader(
    sequence: str,
    psfm_matrix: np.ndarray,
    repo_root: Path,
    dmvfl_root: Path,
    step_timeout_seconds: int,
) -> np.ndarray:
    cache_dir = repo_root / ".analysis-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    output_root = dmvfl_root / "Util" / "database" / "ProtChain"
    with tempfile.TemporaryDirectory(prefix="dmvfl-threader-", dir=cache_dir) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        psfm_path = temp_dir / "input.psfm"
        np.savetxt(psfm_path, psfm_matrix, fmt="%.6f")
        output_base = temp_dir / "thread"

        config_text = (dmvfl_root / "Config.properties").read_text(encoding="utf-8")
        config_text = config_text.replace(
            "/data0/junh/stu/xueqiangf/SPRSA/Util/ProtChain",
            str(output_root.resolve()),
        )
        (temp_dir / "Config.properties").write_text(config_text, encoding="utf-8")

        subprocess.run(
            [
                "java",
                "-jar",
                str((dmvfl_root / "JPSFMThreader.jar").resolve()),
                "query",
                sequence,
                "0.5",
                "1",
                str(psfm_path.resolve()),
                str(output_base.resolve()),
            ],
            cwd=temp_dir,
            check=True,
            capture_output=True,
            text=True,
            timeout=step_timeout_seconds,
        )

        sa_path = output_base.with_suffix(".sa")
        rows = []
        for raw_line in sa_path.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            _position, _residue, _pred, buried, intermediate, exposed = line.split()
            rows.append([float(buried), float(intermediate), float(exposed)])
    return np.array(rows, dtype=np.float32)


def _run_model(
    pssm_matrix: np.ndarray,
    psfm_matrix: np.ndarray,
    pss_matrix: np.ndarray,
    jpsfm_matrix: np.ndarray,
    dmvfl_root: Path,
) -> np.ndarray:
    sys.path.insert(0, str(dmvfl_root.resolve()))
    if "numba" not in sys.modules:
        fake_numba = types.ModuleType("numba")
        fake_numba.jit = lambda func=None, *args, **kwargs: (lambda inner: inner)(func) if callable(func) else (lambda inner: inner)
        sys.modules["numba"] = fake_numba
    if not hasattr(np, "float"):
        np.float = float

    from BiLSTM_SE_Net import LSTMMergeSENet

    model = LSTMMergeSENet()
    state = torch.load(dmvfl_root / "save_model" / "save_model" / "epoch_50", map_location="cpu")
    model.load_state_dict(state)
    model.eval()

    direct = [
        torch.unsqueeze(torch.tensor(_window_matrix(matrix), dtype=torch.float32), 0)
        for matrix in (pssm_matrix, psfm_matrix, pss_matrix, jpsfm_matrix)
    ]
    reverse = [torch.flip(matrix, dims=(2,)) for matrix in direct]

    with torch.no_grad():
        predict00 = model(*direct)
        predict01 = model(*reverse)
        prediction = ((predict00[4] + predict01[4]) / 2.0).squeeze(1).numpy()
    return np.clip(prediction, 0.0, 1.0)


def _window_matrix(matrix: np.ndarray, win_size: int = 15) -> np.ndarray:
    stride = win_size // 2
    padded = np.concatenate([matrix[:stride], matrix, matrix[-stride:]], axis=0)
    windows = np.zeros((matrix.shape[0], win_size * matrix.shape[1]), dtype=np.float32)
    for index in range(stride, matrix.shape[0] + stride):
        windows[index - stride, :] = padded[index - stride : index + stride + 1, :].reshape(-1)
    return windows


def _build_payload(sequence: str, rsa_values: np.ndarray, asa_values: np.ndarray, variant_position: int, args) -> dict:
    exposure_labels = []
    for rsa in rsa_values:
        if rsa >= RSA_EXPOSED_THRESHOLD:
            exposure_labels.append("Exposed")
        elif rsa <= RSA_BURIED_THRESHOLD:
            exposure_labels.append("Buried")
        else:
            exposure_labels.append("Intermediate")

    max_asa = float(np.max(asa_values)) if len(asa_values) else 1.0
    variant_index = variant_position - 1
    variant_rsa = float(rsa_values[variant_index])
    variant_asa = float(asa_values[variant_index])
    variant_label = exposure_labels[variant_index]

    rows = []
    for index, residue in enumerate(sequence, start=1):
        rows.append(
            {
                "Position": index,
                "AA": residue,
                "RSA": round(float(rsa_values[index - 1]), 4),
                "ASA": round(float(asa_values[index - 1]), 2),
                "Exposure": exposure_labels[index - 1],
            }
        )

    interval_features = []
    interval_features.extend(_accessibility_segments(rsa_values, variant_position, threshold=RSA_EXPOSED_THRESHOLD, mode="exposed"))
    interval_features.extend(_accessibility_segments(rsa_values, variant_position, threshold=RSA_BURIED_THRESHOLD, mode="buried"))

    return {
        "summary": f"DMVFL-RSA completed for {len(sequence)} residues; variant-site RSA={variant_rsa:.3f}, ASA={variant_asa:.1f}.",
        "interval_features": interval_features,
        "site_features": [
            site_feature(
                source="DMVFL-RSA",
                label=f"{variant_label} residue",
                category="structure",
                position=variant_position,
                variant_position=variant_position,
                type_name="RSA_ASA_SITE",
                description=(
                    f"Predicted RSA {variant_rsa:.3f} and ASA {variant_asa:.1f} A^2 at residue {variant_position}."
                ),
                evidence=f"Exposure class {variant_label}",
            )
        ],
        "residue_tracks": [
            {
                "key": "dmvfl_rsa",
                "label": "DMVFL-RSA RSA",
                "values": [{"position": index, "value": round(float(value), 4)} for index, value in enumerate(rsa_values, start=1)],
            },
            {
                "key": "dmvfl_asa_scaled",
                "label": "DMVFL-RSA ASA (scaled)",
                "values": [
                    {"position": index, "value": round(float(value) / max_asa, 4)}
                    for index, value in enumerate(asa_values, start=1)
                ],
            },
        ],
        "evidence_rows": [
            {
                "source": "DMVFL-RSA",
                "kind": "rsa_asa",
                "label": "Variant-site accessibility",
                "position": variant_position,
                "summary": f"Predicted {variant_label.lower()} residue with RSA {variant_rsa:.3f} and ASA {variant_asa:.1f} A^2.",
                "score": f"{variant_rsa:.3f}",
            }
        ],
        "tables": [
            table(
                "DMVFL-RSA accessibility summary",
                ["Metric", "Value"],
                [
                    {"Metric": "Mean RSA", "Value": f"{float(np.mean(rsa_values)):.3f}"},
                    {"Metric": "Mean ASA", "Value": f"{float(np.mean(asa_values)):.1f}"},
                    {"Metric": "Exposed residues", "Value": exposure_labels.count("Exposed")},
                    {"Metric": "Buried residues", "Value": exposure_labels.count("Buried")},
                    {"Metric": "Variant-site class", "Value": variant_label},
                ],
            ),
            table("DMVFL-RSA residue accessibility", ["Position", "AA", "RSA", "ASA", "Exposure"], rows),
        ],
    }


def _accessibility_segments(rsa_values: np.ndarray, variant_position: int, *, threshold: float, mode: str) -> list[dict]:
    features: list[dict] = []
    start = None
    for index, value in enumerate(rsa_values, start=1):
        hit = value >= threshold if mode == "exposed" else value <= threshold
        if hit and start is None:
            start = index
        if not hit and start is not None:
            if index - start >= 4:
                features.append(_segment_feature(start, index - 1, mode, variant_position))
            start = None
    if start is not None and len(rsa_values) - start + 1 >= 4:
        features.append(_segment_feature(start, len(rsa_values), mode, variant_position))
    return features


def _segment_feature(start: int, end: int, mode: str, variant_position: int) -> dict:
    if mode == "exposed":
        return interval_feature(
            source="DMVFL-RSA",
            label="Predicted exposed stretch",
            category="structure",
            start=start,
            end=end,
            variant_position=variant_position,
            type_name="RSA_EXPOSED",
            description="Consecutive residues predicted to have high solvent accessibility.",
        )
    return interval_feature(
        source="DMVFL-RSA",
        label="Predicted buried stretch",
        category="structure",
        start=start,
        end=end,
        variant_position=variant_position,
        type_name="RSA_BURIED",
        description="Consecutive residues predicted to have low solvent accessibility.",
    )


if __name__ == "__main__":
    main()
