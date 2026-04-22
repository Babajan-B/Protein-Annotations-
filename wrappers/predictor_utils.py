from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time
from typing import Any
from xml.etree import ElementTree as ET

import requests


DEFAULT_REQUEST_TIMEOUT = 30.0
DEFAULT_POLL_INTERVAL_SECONDS = 2.0
DEFAULT_MAX_WAIT_SECONDS = 180.0
DEFAULT_EMAIL = "academic@example.org"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class PredictorRuntimeError(RuntimeError):
    """Raised when a predictor wrapper fails."""


def build_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--variant-position", required=True, type=int)
    parser.add_argument("--sequence-id", default="")
    parser.add_argument("--gene-symbol", default="")
    parser.add_argument("--wild-type", default="")
    parser.add_argument("--mutant", default="")
    parser.add_argument("--request-timeout", type=float, default=DEFAULT_REQUEST_TIMEOUT)
    parser.add_argument("--max-wait-seconds", type=float, default=DEFAULT_MAX_WAIT_SECONDS)
    parser.add_argument("--poll-interval-seconds", type=float, default=DEFAULT_POLL_INTERVAL_SECONDS)
    return parser


def read_single_fasta(path: str) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            continue
        sequence_parts.append(line)
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise PredictorRuntimeError(f"No protein sequence found in {path}.")
    return header, sequence


def write_payload(path: str, payload: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def interval_feature(
    *,
    source: str,
    label: str,
    category: str,
    start: int,
    end: int,
    variant_position: int,
    type_name: str = "",
    description: str = "",
    evidence: str = "",
) -> dict[str, Any]:
    return {
        "source": source,
        "label": label,
        "type": type_name or label,
        "category": category,
        "start": start,
        "end": end,
        "description": description or label,
        "evidence": evidence,
        "variant_hit": start <= variant_position <= end,
    }


def site_feature(
    *,
    source: str,
    label: str,
    category: str,
    position: int,
    variant_position: int,
    type_name: str = "",
    description: str = "",
    evidence: str = "",
) -> dict[str, Any]:
    return {
        "source": source,
        "label": label,
        "type": type_name or label,
        "category": category,
        "position": position,
        "description": description or label,
        "evidence": evidence,
        "variant_hit": position == variant_position,
    }


def table(title: str, columns: list[str], rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {"title": title, "columns": columns, "rows": rows}


class EbiJobClient:
    def __init__(
        self,
        tool: str,
        *,
        request_timeout: float,
        max_wait_seconds: float,
        poll_interval_seconds: float,
        email: str | None = None,
    ) -> None:
        self.tool = tool
        self.base_url = f"https://www.ebi.ac.uk/Tools/services/rest/{tool}"
        self.request_timeout = request_timeout
        self.max_wait_seconds = max_wait_seconds
        self.poll_interval_seconds = poll_interval_seconds
        self.email = email or os.environ.get("EBI_JOB_EMAIL", DEFAULT_EMAIL)
        self.headers = {"User-Agent": "ProteinMutationWorkbench/0.2 academic prototype"}

    def submit(self, params: dict[str, Any]) -> str:
        payload = {"email": self.email, **params}
        response = requests.post(
            f"{self.base_url}/run",
            data=payload,
            headers=self.headers,
            timeout=self.request_timeout,
        )
        if response.status_code != 200:
            raise PredictorRuntimeError(
                f"{self.tool} submission failed with HTTP {response.status_code}: {response.text[:300].strip()}"
            )
        return response.text.strip()

    def wait(self, job_id: str) -> None:
        deadline = time.monotonic() + self.max_wait_seconds
        while time.monotonic() < deadline:
            response = requests.get(
                f"{self.base_url}/status/{job_id}",
                headers=self.headers,
                timeout=self.request_timeout,
            )
            if response.status_code != 200:
                raise PredictorRuntimeError(
                    f"{self.tool} status check failed with HTTP {response.status_code}: {response.text[:200].strip()}"
                )
            status = response.text.strip().upper()
            if status == "FINISHED":
                return
            if status in {"ERROR", "FAILURE", "FAILED", "NOT_FOUND"}:
                raise PredictorRuntimeError(f"{self.tool} job {job_id} ended with status {status}.")
            time.sleep(self.poll_interval_seconds)
        raise PredictorRuntimeError(f"{self.tool} job {job_id} did not finish within {self.max_wait_seconds:.0f}s.")

    def result_types(self, job_id: str) -> list[str]:
        response = requests.get(
            f"{self.base_url}/resulttypes/{job_id}",
            headers=self.headers,
            timeout=self.request_timeout,
        )
        if response.status_code != 200:
            raise PredictorRuntimeError(
                f"{self.tool} result-type query failed with HTTP {response.status_code}: {response.text[:200].strip()}"
            )
        root = ET.fromstring(response.text)
        return [node.text for node in root.findall(".//identifier") if node.text]

    def result(self, job_id: str, result_type: str) -> str:
        response = requests.get(
            f"{self.base_url}/result/{job_id}/{result_type}",
            headers=self.headers,
            timeout=self.request_timeout,
        )
        if response.status_code != 200:
            raise PredictorRuntimeError(
                f"{self.tool} result fetch for {result_type} failed with HTTP {response.status_code}: "
                f"{response.text[:200].strip()}"
            )
        return response.text

    def default_notices(self) -> list[str]:
        if self.email == DEFAULT_EMAIL:
            return [
                "EBI job wrappers are using the placeholder email academic@example.org. Set EBI_JOB_EMAIL for production use."
            ]
        return []
