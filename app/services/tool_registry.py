from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shlex
import sys


@dataclass(frozen=True)
class PredictorSpec:
    key: str
    label: str
    category: str
    summary: str
    command_env: str
    docs_url: str
    resource_hint: str


ROOT_DIR = Path(__file__).resolve().parents[2]
WRAPPERS_DIR = ROOT_DIR / "wrappers"
BIO_TOOLS_DIR = ROOT_DIR / "envs" / "bio-tools"


PREDICTOR_SPECS = [
    PredictorSpec(
        key="interproscan",
        label="InterProScan",
        category="domains",
        summary="Sequence domain and motif mapping across InterPro member databases.",
        command_env="INTERPROSCAN_COMMAND_TEMPLATE",
        docs_url="https://interpro-documentation.readthedocs.io/en/latest/download.html",
        resource_hint="Defaults to the official EBI iprscan5 web-service wrapper unless you override it with a local InterProScan install.",
    ),
    PredictorSpec(
        key="dmvfl_rsa",
        label="DMVFL-RSA",
        category="structure",
        summary="Relative and absolute solvent accessibility prediction using the bundled DMVFL-RSA model and local ProtChain resources.",
        command_env="DMVFL_RSA_COMMAND_TEMPLATE",
        docs_url="https://github.com/XueqiangF/DMVFL-RSA",
        resource_hint="Defaults to the local DMVFL-RSA-main repository in this workspace when its model assets are present.",
    ),
    PredictorSpec(
        key="hmmer",
        label="HMMER phmmer",
        category="conservation",
        summary="Sequence-homology search against SwissProt plus conservation scoring from the returned alignment.",
        command_env="HMMER_COMMAND_TEMPLATE",
        docs_url="http://hmmer.org/documentation.html",
        resource_hint="Defaults to the official EBI HMMER phmmer service unless you override it with a local database-backed workflow.",
    ),
    PredictorSpec(
        key="hhblits",
        label="HHblits",
        category="conservation",
        summary="Fast profile-profile search for deep alignments and conservation context.",
        command_env="HHBLITS_COMMAND_TEMPLATE",
        docs_url="https://github.com/soedinglab/hh-suite",
        resource_hint="Requires HH-suite plus a configured Uniclust or UniRef-style database in the wrapper command.",
    ),
    PredictorSpec(
        key="psipred",
        label="PSIPRED",
        category="structure",
        summary="Secondary structure prediction from sequence, using PSIPRED single-sequence mode by default.",
        command_env="PSIPRED_COMMAND_TEMPLATE",
        docs_url="https://github.com/psipred/psipred",
        resource_hint="Defaults to the local PSIPRED binaries in envs/bio-tools when present.",
    ),
    PredictorSpec(
        key="phobius",
        label="Phobius",
        category="topology",
        summary="Signal peptide and membrane topology prediction from sequence.",
        command_env="PHOBIUS_COMMAND_TEMPLATE",
        docs_url="https://www.ebi.ac.uk/jdispatcher/seqstats/phobius",
        resource_hint="Defaults to the official EBI Phobius web-service wrapper unless you override it.",
    ),
    PredictorSpec(
        key="deeptmhmm",
        label="DeepTMHMM",
        category="topology",
        summary="Transmembrane helix and topology prediction.",
        command_env="DEEPTMHMM_COMMAND_TEMPLATE",
        docs_url="https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0",
        resource_hint="Use a local package or wrapper around the official standalone resources if available to you.",
    ),
    PredictorSpec(
        key="signalp6",
        label="SignalP 6.0",
        category="topology",
        summary="Signal peptide detection and cleavage-site prediction.",
        command_env="SIGNALP6_COMMAND_TEMPLATE",
        docs_url="https://services.healthtech.dtu.dk/services/SignalP-6.0/",
        resource_hint="Requires the academic download and model resources, surfaced via a wrapper command.",
    ),
    PredictorSpec(
        key="deeploc",
        label="DeepLoc 2.x",
        category="localization",
        summary="Subcellular localization and sorting-signal prediction for eukaryotic proteins.",
        command_env="DEEPLOC_COMMAND_TEMPLATE",
        docs_url="https://services.healthtech.dtu.dk/services/DeepLoc-2.0/",
        resource_hint="Requires the DTU package download or an equivalent wrapper command on a prepared environment.",
    ),
    PredictorSpec(
        key="iupred3",
        label="IUPred3",
        category="disorder",
        summary="Intrinsic disorder and binding-region prediction from sequence.",
        command_env="IUPRED3_COMMAND_TEMPLATE",
        docs_url="https://iupred3.elte.hu/download_new",
        resource_hint="Academic download only; provide the local command via the environment once received.",
    ),
    PredictorSpec(
        key="disopred3",
        label="DISOPRED3",
        category="disorder",
        summary="Intrinsic disorder and disordered-binding-site prediction.",
        command_env="DISOPRED3_COMMAND_TEMPLATE",
        docs_url="https://github.com/psipred/disopred",
        resource_hint="Requires a DISOPRED installation and its supporting profile pipeline.",
    ),
    PredictorSpec(
        key="musitedeep",
        label="MusiteDeep",
        category="ptm",
        summary="Sequence-based PTM prediction, especially human phosphorylation models.",
        command_env="MUSITEDEEP_COMMAND_TEMPLATE",
        docs_url="https://github.com/duolinwang/MusiteDeep",
        resource_hint="Requires the MusiteDeep environment and model weights exposed through a wrapper command.",
    ),
    PredictorSpec(
        key="evcouplings",
        label="EVcouplings",
        category="conservation",
        summary="Co-evolutionary coupling analysis for residue interaction context.",
        command_env="EVCOUPLINGS_COMMAND_TEMPLATE",
        docs_url="https://evcouplings.readthedocs.io/",
        resource_hint="Requires EVcouplings plus alignment resources and often additional databases or plmc binaries.",
    ),
    PredictorSpec(
        key="consurf",
        label="ConSurf",
        category="conservation",
        summary="Residue conservation scoring from homologous sequence alignments.",
        command_env="CONSURF_COMMAND_TEMPLATE",
        docs_url="https://consurf.tau.ac.il/",
        resource_hint="Provide a local wrapper if you have an academic installation or an internal equivalent conservation workflow.",
    ),
    PredictorSpec(
        key="netsurfp3",
        label="NetSurfP-3.0",
        category="structure",
        summary="Per-residue secondary structure (H/E/C), relative solvent accessibility, backbone torsion angles, and disorder probability predicted from sequence using deep learning.",
        command_env="NETSURFP3_COMMAND_TEMPLATE",
        docs_url="https://services.healthtech.dtu.dk/services/NetSurfP-3.0/",
        resource_hint="Requires a local NetSurfP-3.0 installation. Download from DTU Health Tech, then configure via NETSURFP3_COMMAND_TEMPLATE or place the netsurfp3 binary in envs/bio-tools/bin.",
    ),
]


def resolve_command_template(spec: PredictorSpec) -> str:
    explicit = os.environ.get(spec.command_env, "").strip()
    if explicit:
        return explicit
    return _default_command_template(spec.key)


def build_resource_cards(enabled_keys: set[str] | None = None) -> list[dict[str, str]]:
    cards = []
    for spec in PREDICTOR_SPECS:
        if enabled_keys and spec.key not in enabled_keys:
            continue
        configured = bool(resolve_command_template(spec))
        cards.append(
            {
                "key": spec.key,
                "label": spec.label,
                "category": spec.category,
                "summary": spec.summary,
                "command_env": spec.command_env,
                "configured": "configured" if configured else "missing",
                "docs_url": spec.docs_url,
                "resource_hint": spec.resource_hint,
            }
        )
    return cards


def _default_command_template(key: str) -> str:
    python_bin = _default_python_executable()
    if not python_bin:
        return ""

    commands = {
        "interproscan": _wrapper_command(
            python_bin,
            "run_iprscan5_remote.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position} "
            "--sequence-id {sequence_id} --max-wait-seconds 180",
        ),
        "dmvfl_rsa": _wrapper_command(
            python_bin,
            "run_dmvfl_rsa.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position} "
            "--sequence-id {sequence_id} --max-wait-seconds 180",
            requires=[
                ROOT_DIR / "DMVFL-RSA-main" / "save_model" / "save_model" / "epoch_50",
                ROOT_DIR / "DMVFL-RSA-main" / "Util" / "database" / "ProtChain",
                BIO_TOOLS_DIR / "bin" / "psipred",
                BIO_TOOLS_DIR / "bin" / "seq2mtx",
                BIO_TOOLS_DIR / "bin" / "psipass2",
            ],
        ),
        "hmmer": _wrapper_command(
            python_bin,
            "run_hmmer_phmmer_remote.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position} "
            "--max-wait-seconds 120 --nhits 64",
        ),
        "psipred": _wrapper_command(
            python_bin,
            "run_psipred_single.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position}",
            requires=[BIO_TOOLS_DIR / "bin" / "psipred", BIO_TOOLS_DIR / "bin" / "seq2mtx", BIO_TOOLS_DIR / "bin" / "psipass2"],
        ),
        "phobius": _wrapper_command(
            python_bin,
            "run_phobius_remote.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position} "
            "--max-wait-seconds 120",
        ),
        "netsurfp3": _wrapper_command(
            python_bin,
            "run_netsurfp3_remote.py",
            "--input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position}",
            requires=[BIO_TOOLS_DIR / "bin" / "netsurfp3"],
        ),
    }
    return commands.get(key, "")


def _wrapper_command(python_bin: str, wrapper_name: str, args: str, *, requires: list[Path] | None = None) -> str:
    wrapper_path = WRAPPERS_DIR / wrapper_name
    if not wrapper_path.exists():
        return ""
    if requires and any(not path.exists() for path in requires):
        return ""
    return f"{shlex.quote(python_bin)} {shlex.quote(str(wrapper_path))} {args}"


def _default_python_executable() -> str:
    candidates = [
        ROOT_DIR / ".venv" / "bin" / "python",
        ROOT_DIR / "envs" / "bio-tools" / "bin" / "python",
        Path(sys.executable),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return ""
