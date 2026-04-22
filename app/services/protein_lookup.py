from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any

import requests


ACCESSION_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9](?:-\d+)?$|^[A-NR-Z][0-9]{5}(?:-\d+)?$", re.IGNORECASE)


class ProteinLookupError(RuntimeError):
    """Raised when a human protein cannot be resolved from the public data source."""


@dataclass(frozen=True)
class ProteinRecord:
    accession: str
    entry_name: str
    gene_symbol: str
    protein_name: str
    organism_name: str
    sequence: str
    features: list[dict[str, Any]]
    refseq_proteins: list[str]
    synonyms: list[str]

    @property
    def length(self) -> int:
        return len(self.sequence)


class HumanProteinResolver:
    def __init__(self, base_url: str, timeout_seconds: float) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.session = requests.Session()

    def resolve(self, query: str, isoform: str | None = None) -> ProteinRecord:
        normalized_query = (isoform or query or "").strip()
        if not normalized_query:
            raise ProteinLookupError("A gene symbol, protein name, or UniProt accession is required.")

        accession = isoform.strip() if isoform else self._search_accession(normalized_query)
        payload = self._request_json(f"{self.base_url}/uniprotkb/{accession}.json")

        try:
            gene_symbol = payload["genes"][0]["geneName"]["value"]
        except (KeyError, IndexError):
            gene_symbol = ""

        protein_name = self._extract_protein_name(payload)
        sequence = payload["sequence"]["value"]
        organism_name = payload["organism"]["scientificName"]
        features = payload.get("features", [])
        refseq_proteins = self._extract_refseq_proteins(payload)
        synonyms = self._extract_synonyms(payload)

        return ProteinRecord(
            accession=payload["primaryAccession"],
            entry_name=payload["uniProtkbId"],
            gene_symbol=gene_symbol,
            protein_name=protein_name,
            organism_name=organism_name,
            sequence=sequence,
            features=features,
            refseq_proteins=refseq_proteins,
            synonyms=synonyms,
        )

    def _search_accession(self, query: str) -> str:
        if ACCESSION_RE.match(query):
            return query.upper()

        escaped_query = query.replace('"', "").strip()
        quoted_query = f'"{escaped_query}"'
        search_query = (
            f"((gene_exact:{quoted_query}) OR (gene:{quoted_query}) OR (protein_name:{quoted_query})) "
            "AND organism_id:9606"
        )
        params = {
            "query": search_query,
            "size": 5,
            "fields": "accession,id,protein_name,gene_names,organism_name",
            "format": "json",
        }
        payload = self._request_json(f"{self.base_url}/uniprotkb/search", params=params)
        results = payload.get("results", [])
        if not results:
            raise ProteinLookupError(f"No human UniProt entry was found for '{query}'.")

        first = results[0]
        accession = first.get("primaryAccession")
        if not accession:
            raise ProteinLookupError(f"UniProt search returned an entry without an accession for '{query}'.")
        return accession

    def _request_json(self, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        try:
            response = self.session.get(url, params=params, timeout=self.timeout_seconds)
            response.raise_for_status()
        except requests.RequestException as exc:
            raise ProteinLookupError(f"Protein lookup failed: {exc}") from exc
        return response.json()

    @staticmethod
    def _extract_protein_name(payload: dict[str, Any]) -> str:
        description = payload.get("proteinDescription", {})

        recommended = description.get("recommendedName", {})
        full_name = recommended.get("fullName", {}).get("value")
        if full_name:
            return full_name

        submission_names = description.get("submissionNames", [])
        if submission_names:
            submitted = submission_names[0].get("fullName", {}).get("value")
            if submitted:
                return submitted

        return payload.get("uniProtkbId", "Unknown protein")

    @staticmethod
    def _extract_refseq_proteins(payload: dict[str, Any]) -> list[str]:
        proteins = []
        for reference in payload.get("uniProtKBCrossReferences", []):
            if reference.get("database") != "RefSeq":
                continue
            accession = reference.get("id")
            if accession:
                proteins.append(accession.split(".")[0])
        return proteins

    @staticmethod
    def _extract_synonyms(payload: dict[str, Any]) -> list[str]:
        description = payload.get("proteinDescription", {})
        recommended = description.get("recommendedName", {})
        names = []

        recommended_full = recommended.get("fullName", {}).get("value")
        if recommended_full:
            names.append(recommended_full)

        for alternative in description.get("alternativeNames", []):
            value = alternative.get("fullName", {}).get("value")
            if value:
                names.append(value)

        return names
