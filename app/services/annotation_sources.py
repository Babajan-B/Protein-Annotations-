from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import requests

from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import MissenseVariant


VISIBLE_FEATURE_TYPES = {
    "DOMAIN",
    "REGION",
    "REPEAT",
    "MOTIF",
    "ZN_FING",
    "DNA_BIND",
    "ACT_SITE",
    "BINDING",
    "SITE",
    "METAL",
    "MOD_RES",
    "LIPID",
    "CARBOHYD",
    "DISULFID",
    "SIGNAL",
    "TRANSIT",
    "TRANSMEM",
    "INTRAMEM",
    "TOPO_DOM",
    "HELIX",
    "STRAND",
    "TURN",
    "COILED",
    "COMPBIAS",
    "VARIANT",
    "MUTAGEN",
}

TRACK_CATEGORY_MAP = {
    "DOMAIN": "domains",
    "REGION": "domains",
    "REPEAT": "domains",
    "MOTIF": "domains",
    "ZN_FING": "domains",
    "DNA_BIND": "domains",
    "SIGNAL": "topology",
    "TRANSIT": "topology",
    "TRANSMEM": "topology",
    "INTRAMEM": "topology",
    "TOPO_DOM": "topology",
    "COILED": "structure",
    "COMPBIAS": "structure",
    "HELIX": "structure",
    "STRAND": "structure",
    "TURN": "structure",
    "ACT_SITE": "sites",
    "BINDING": "sites",
    "SITE": "sites",
    "METAL": "sites",
    "MOD_RES": "ptm",
    "LIPID": "ptm",
    "CARBOHYD": "ptm",
    "DISULFID": "ptm",
    "VARIANT": "evidence",
    "MUTAGEN": "evidence",
}

CATEGORY_LABELS = {
    "domains": "Domains and motifs",
    "topology": "Topology and targeting",
    "structure": "Curated structure features",
    "sites": "Functional sites",
    "ptm": "PTMs and covalent features",
    "evidence": "Known variants and mutagenesis",
}

CATEGORY_COLORS = {
    "domains": "#2c8f71",
    "topology": "#4b6cb7",
    "structure": "#9a6df0",
    "sites": "#d86f45",
    "ptm": "#c08a1f",
    "evidence": "#b53f66",
}


@dataclass(frozen=True)
class SourceStatus:
    source: str
    status: str
    summary: str


class AnnotationAggregator:
    def __init__(
        self,
        *,
        interpro_base_url: str,
        ebi_proteins_base_url: str,
        ensembl_base_url: str,
        ncbi_base_url: str,
        mavedb_base_url: str,
        timeout_seconds: float,
        user_agent: str,
    ) -> None:
        self.interpro_base_url = interpro_base_url.rstrip("/")
        self.ebi_proteins_base_url = ebi_proteins_base_url.rstrip("/")
        self.ensembl_base_url = ensembl_base_url.rstrip("/")
        self.ncbi_base_url = ncbi_base_url.rstrip("/")
        self.mavedb_base_url = mavedb_base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.session = requests.Session()
        self.default_headers = {
            "Accept": "application/json",
            "User-Agent": user_agent,
        }

    def collect(self, protein: ProteinRecord, variant: MissenseVariant) -> dict[str, Any]:
        interval_features, site_features, evidence_rows = self._collect_uniprot_features(protein, variant)

        source_statuses = [
            SourceStatus(
                source="UniProt curated features",
                status="ok",
                summary=f"Loaded {len(interval_features)} interval features and {len(site_features)} residue-linked features.",
            )
        ]
        notices: list[str] = []

        interpro_features, interpro_status = self._collect_interpro_features(protein.accession)
        if interpro_features:
            interval_features.extend(interpro_features)
        source_statuses.append(interpro_status)
        if interpro_status.status != "ok":
            notices.append(interpro_status.summary)

        ebi_rows, ebi_status = self._collect_ebi_proteins_features(protein.accession, variant.position)
        if ebi_rows:
            evidence_rows.extend(ebi_rows)
        source_statuses.append(ebi_status)
        if ebi_status.status != "ok":
            notices.append(ebi_status.summary)

        vep_rows, vep_status = self._collect_vep_rows(protein, variant)
        if vep_rows:
            evidence_rows.extend(vep_rows)
        source_statuses.append(vep_status)
        if vep_status.status != "ok":
            notices.append(vep_status.summary)

        clinvar_rows, clinvar_status = self._collect_clinvar_rows(protein, variant)
        if clinvar_rows:
            evidence_rows.extend(clinvar_rows)
        source_statuses.append(clinvar_status)
        if clinvar_status.status != "ok":
            notices.append(clinvar_status.summary)

        mavedb_rows, mavedb_status = self._collect_mavedb_rows(protein, variant)
        if mavedb_rows:
            evidence_rows.extend(mavedb_rows)
        source_statuses.append(mavedb_status)
        if mavedb_status.status != "ok":
            notices.append(mavedb_status.summary)

        architecture_rows = self._build_architecture_rows(interval_features, protein.length, variant.position)
        feature_counts = self._build_feature_counts(interval_features, site_features)
        evidence_plot = self._build_evidence_plot(evidence_rows, protein.length, variant.position)

        return {
            "interval_features": sorted(interval_features, key=lambda row: (row["start"], row["end"], row["source"])),
            "site_features": sorted(site_features, key=lambda row: (row["position"], row["source"], row["label"])),
            "evidence_rows": evidence_rows,
            "architecture_rows": architecture_rows,
            "feature_counts": feature_counts,
            "evidence_plot": evidence_plot,
            "source_status": [status.__dict__ for status in source_statuses],
            "notices": sorted(set(notices)),
        }

    def _collect_uniprot_features(
        self, protein: ProteinRecord, variant: MissenseVariant
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
        interval_features: list[dict[str, Any]] = []
        site_features: list[dict[str, Any]] = []
        evidence_rows: list[dict[str, Any]] = []

        for feature in protein.features:
            feature_type = feature.get("type", "").upper()
            if feature_type not in VISIBLE_FEATURE_TYPES:
                continue

            location = feature.get("location", {})
            start = _extract_position(location.get("start")) or _extract_position(location.get("begin"))
            end = _extract_position(location.get("end")) or _extract_position(location.get("stop"))
            if not start and not end:
                continue

            if start is None:
                start = end
            if end is None:
                end = start

            row = {
                "source": "UniProt",
                "type": feature_type,
                "category": TRACK_CATEGORY_MAP.get(feature_type, "domains"),
                "label": feature.get("description") or feature_type.replace("_", " ").title(),
                "description": feature.get("description") or feature_type.replace("_", " ").title(),
                "start": start,
                "end": end,
                "evidence": _format_evidence_codes(feature.get("evidences", [])),
                "variant_hit": start <= variant.position <= end,
            }

            if start == end:
                row["position"] = start
                site_features.append(row)
            else:
                interval_features.append(row)

            if feature_type in {"VARIANT", "MUTAGEN"}:
                evidence_rows.append(
                    {
                        "source": "UniProt",
                        "kind": "curated_variant",
                        "label": row["label"],
                        "position": start,
                        "summary": row["description"],
                        "link_hint": protein.accession,
                        "score": "",
                    }
                )

        return interval_features, site_features, evidence_rows

    def _collect_interpro_features(self, accession: str) -> tuple[list[dict[str, Any]], SourceStatus]:
        candidate_urls = [
            f"{self.interpro_base_url}/entry/interpro/protein/uniprot/{accession}",
            f"{self.interpro_base_url}/entry/all/protein/uniprot/{accession}",
        ]
        for url in candidate_urls:
            payload = self._request_json(url)
            if not payload:
                continue

            results = payload.get("results") if isinstance(payload, dict) else payload
            if not isinstance(results, list):
                continue

            features: list[dict[str, Any]] = []
            for result in results:
                metadata = result.get("metadata", {}) if isinstance(result, dict) else {}
                label = metadata.get("name") or metadata.get("accession") or "InterPro feature"
                source_database = metadata.get("source_database") or "InterPro"

                for fragment in _walk_interpro_fragments(result):
                    start = fragment.get("start")
                    end = fragment.get("end")
                    if not start or not end:
                        continue
                    features.append(
                        {
                            "source": f"InterPro:{source_database}",
                            "type": "INTERPRO",
                            "category": "domains",
                            "label": label,
                            "description": metadata.get("description") or label,
                            "start": start,
                            "end": end,
                            "evidence": metadata.get("accession", ""),
                            "variant_hit": False,
                        }
                    )

            if features:
                return features, SourceStatus(
                    source="InterPro API",
                    status="ok",
                    summary=f"Loaded {len(features)} InterPro intervals.",
                )

        return [], SourceStatus(
            source="InterPro API",
            status="unavailable",
            summary="InterPro enrichment was attempted, but no compatible public API response was available for this entry.",
        )

    def _collect_ebi_proteins_features(
        self, accession: str, variant_position: int
    ) -> tuple[list[dict[str, Any]], SourceStatus]:
        candidate_urls = [
            f"{self.ebi_proteins_base_url}/features/{accession}",
            f"{self.ebi_proteins_base_url}/variation/{accession}",
        ]

        rows: list[dict[str, Any]] = []
        for url in candidate_urls:
            payload = self._request_json(url)
            if not payload:
                continue

            entries = []
            if isinstance(payload, dict):
                entries = payload.get("features") or payload.get("sequenceFeatures") or payload.get("featuresData") or []
            elif isinstance(payload, list):
                entries = payload

            for entry in entries:
                start = _extract_position(entry.get("begin")) or _extract_position(entry.get("start"))
                end = _extract_position(entry.get("end")) or start
                if not start:
                    continue
                if not (start <= variant_position <= (end or start)):
                    continue
                rows.append(
                    {
                        "source": "EBI Proteins API",
                        "kind": "protein_feature",
                        "label": entry.get("type") or entry.get("category") or "Feature",
                        "position": start,
                        "summary": entry.get("description") or entry.get("consequenceType") or "Residue-linked feature",
                        "link_hint": accession,
                        "score": "",
                    }
                )

            if rows:
                return rows, SourceStatus(
                    source="EBI Proteins API",
                    status="ok",
                    summary=f"Loaded {len(rows)} residue-linked EBI proteins annotations.",
                )

        return [], SourceStatus(
            source="EBI Proteins API",
            status="unavailable",
            summary="EBI Proteins enrichment was attempted, but no compatible public feature payload was returned.",
        )

    def _collect_vep_rows(
        self, protein: ProteinRecord, variant: MissenseVariant
    ) -> tuple[list[dict[str, Any]], SourceStatus]:
        accessions = protein.refseq_proteins or []
        if not accessions:
            return [], SourceStatus(
                source="Ensembl VEP",
                status="unavailable",
                summary="No RefSeq protein accession was available from UniProt for VEP protein HGVS submission.",
            )

        rows: list[dict[str, Any]] = []
        for accession in accessions[:3]:
            hgvs = f"{accession}:{variant.hgvs_protein}"
            payload = self._request_json(
                f"{self.ensembl_base_url}/vep/human/hgvs/{hgvs}",
                headers={"Content-Type": "application/json", "Accept": "application/json"},
            )
            if not isinstance(payload, list) or not payload:
                continue

            first = payload[0]
            consequences = first.get("most_severe_consequence") or ""
            colocated = first.get("colocated_variants") or []
            identifier = first.get("id") or (colocated[0].get("id") if colocated else accession)
            rows.append(
                {
                    "source": "Ensembl VEP",
                    "kind": "vep_consequence",
                    "label": identifier,
                    "position": variant.position,
                    "summary": consequences.replace("_", " ") if consequences else "Consequence available",
                    "link_hint": accession,
                    "score": "",
                }
            )
            if colocated:
                for colocated_variant in colocated[:3]:
                    rows.append(
                        {
                            "source": "Ensembl VEP",
                            "kind": "vep_colocated",
                            "label": colocated_variant.get("id") or "Colocated variant",
                            "position": variant.position,
                            "summary": colocated_variant.get("clin_sig", [""])[0] or "Known colocated variant",
                            "link_hint": accession,
                            "score": "",
                        }
                    )

            return rows, SourceStatus(
                source="Ensembl VEP",
                status="ok",
                summary=f"Loaded VEP consequence data using {accession}.",
            )

        return [], SourceStatus(
            source="Ensembl VEP",
            status="unavailable",
            summary="VEP lookup was attempted, but the public endpoint did not return protein HGVS consequence data for the current mapping.",
        )

    def _collect_clinvar_rows(
        self, protein: ProteinRecord, variant: MissenseVariant
    ) -> tuple[list[dict[str, Any]], SourceStatus]:
        if not protein.gene_symbol:
            return [], SourceStatus(
                source="ClinVar",
                status="unavailable",
                summary="ClinVar lookup was skipped because no gene symbol was available on the resolved protein.",
            )

        term = f'{protein.gene_symbol}[gene] AND "{variant.hgvs_protein}"'
        search_payload = self._request_json(
            f"{self.ncbi_base_url}/esearch.fcgi",
            params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 5},
        )
        if not isinstance(search_payload, dict):
            return [], SourceStatus(
                source="ClinVar",
                status="unavailable",
                summary="ClinVar esearch did not return a JSON payload.",
            )

        ids = search_payload.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return [], SourceStatus(
                source="ClinVar",
                status="ok",
                summary="ClinVar lookup completed, but no matching records were found for the submitted protein change.",
            )

        summary_payload = self._request_json(
            f"{self.ncbi_base_url}/esummary.fcgi",
            params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
        )
        result_block = summary_payload.get("result", {}) if isinstance(summary_payload, dict) else {}
        rows: list[dict[str, Any]] = []

        for record_id in ids:
            row = result_block.get(record_id, {})
            title = row.get("title") or row.get("accession") or record_id
            clinical = row.get("clinical_significance", {})
            significance = clinical.get("description") if isinstance(clinical, dict) else ""
            rows.append(
                {
                    "source": "ClinVar",
                    "kind": "clinvar",
                    "label": title,
                    "position": variant.position,
                    "summary": significance or "ClinVar record",
                    "link_hint": row.get("accession") or record_id,
                    "score": "",
                }
            )

        return rows, SourceStatus(
            source="ClinVar",
            status="ok",
            summary=f"Loaded {len(rows)} ClinVar record summaries.",
        )

    def _collect_mavedb_rows(
        self, protein: ProteinRecord, variant: MissenseVariant
    ) -> tuple[list[dict[str, Any]], SourceStatus]:
        query = protein.gene_symbol or protein.accession
        payload = self._request_json(
            f"{self.mavedb_base_url}/score-sets",
            params={"search": query},
        )
        if not payload:
            return [], SourceStatus(
                source="MaveDB",
                status="unavailable",
                summary="MaveDB search did not return a compatible public JSON response for this query.",
            )

        items = payload.get("results") if isinstance(payload, dict) else payload
        if not isinstance(items, list) or not items:
            return [], SourceStatus(
                source="MaveDB",
                status="ok",
                summary="MaveDB search completed, but no score sets were found for the current target.",
            )

        rows = []
        variant_tokens = {variant.short_notation.upper(), variant.hgvs_protein.upper()}
        for item in items[:5]:
            title = item.get("title") or item.get("urn") or "MaveDB score set"
            text = " ".join(str(item.get(key, "")) for key in ("title", "description", "urn")).upper()
            if protein.gene_symbol and protein.gene_symbol.upper() not in text and not any(token in text for token in variant_tokens):
                continue
            rows.append(
                {
                    "source": "MaveDB",
                    "kind": "functional_assay",
                    "label": title,
                    "position": variant.position,
                    "summary": item.get("urn") or "Potential functional assay dataset",
                    "link_hint": item.get("urn") or title,
                    "score": item.get("numVariants", ""),
                }
            )

        return rows, SourceStatus(
            source="MaveDB",
            status="ok" if rows else "ok",
            summary=f"Loaded {len(rows)} MaveDB dataset candidates." if rows else "MaveDB search completed with no obvious variant-specific dataset match.",
        )

    def _build_architecture_rows(
        self, interval_features: list[dict[str, Any]], sequence_length: int, variant_position: int
    ) -> list[dict[str, Any]]:
        grouped: dict[str, list[dict[str, Any]]] = {}
        for feature in interval_features:
            category = feature["category"]
            grouped.setdefault(category, []).append(feature)

        rows = []
        for category in ("domains", "topology", "structure", "sites", "ptm", "evidence"):
            items = grouped.get(category, [])
            rendered_segments = []
            for feature in items[:10]:
                rendered_segments.append(
                    {
                        **feature,
                        "left_pct": round(((feature["start"] - 1) / sequence_length) * 100, 3),
                        "width_pct": max(round(((feature["end"] - feature["start"] + 1) / sequence_length) * 100, 3), 0.7),
                        "color": CATEGORY_COLORS[category],
                    }
                )
            rows.append(
                {
                    "key": category,
                    "label": CATEGORY_LABELS[category],
                    "segments": rendered_segments,
                    "variant_pct": round((variant_position / sequence_length) * 100, 3),
                    "has_content": bool(rendered_segments),
                }
            )
        return rows

    def _build_feature_counts(
        self, interval_features: list[dict[str, Any]], site_features: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        counts: dict[str, int] = {}
        for feature in interval_features + site_features:
            label = CATEGORY_LABELS.get(feature["category"], feature["category"])
            counts[label] = counts.get(label, 0) + 1
        return [{"label": label, "count": count} for label, count in sorted(counts.items())]

    @staticmethod
    def _build_evidence_plot(
        evidence_rows: list[dict[str, Any]], sequence_length: int, variant_position: int
    ) -> list[dict[str, Any]]:
        points = []
        for row in evidence_rows:
            position = row.get("position")
            if not position:
                continue
            points.append(
                {
                    **row,
                    "left_pct": round((position / sequence_length) * 100, 3),
                    "hit_variant_site": position == variant_position,
                }
            )
        return points

    def _request_json(
        self,
        url: str,
        params: dict[str, Any] | None = None,
        headers: dict[str, str] | None = None,
    ) -> dict[str, Any] | list[Any] | None:
        merged_headers = dict(self.default_headers)
        if headers:
            merged_headers.update(headers)

        try:
            response = self.session.get(
                url,
                params=params,
                headers=merged_headers,
                timeout=self.timeout_seconds,
            )
            response.raise_for_status()
        except requests.RequestException:
            return None

        try:
            return response.json()
        except ValueError:
            return None


def _extract_position(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    if isinstance(value, dict):
        raw = value.get("value") or value.get("position") or value.get("start") or value.get("end")
        if isinstance(raw, int):
            return raw
        if isinstance(raw, str) and raw.isdigit():
            return int(raw)
    return None


def _format_evidence_codes(entries: list[dict[str, Any]]) -> str:
    codes = []
    for entry in entries:
        code = entry.get("evidenceCode")
        if code:
            codes.append(code)
    return ", ".join(codes)


def _walk_interpro_fragments(result: dict[str, Any]) -> list[dict[str, int]]:
    fragments: list[dict[str, int]] = []
    for protein_entry in result.get("proteins", []):
        for location in protein_entry.get("entry_protein_locations", []) or protein_entry.get("locations", []):
            for fragment in location.get("fragments", []):
                start = _extract_position(fragment.get("start"))
                end = _extract_position(fragment.get("end"))
                if start and end:
                    fragments.append({"start": start, "end": end})

    if fragments:
        return fragments

    for location in result.get("entry_protein_locations", []):
        for fragment in location.get("fragments", []):
            start = _extract_position(fragment.get("start"))
            end = _extract_position(fragment.get("end"))
            if start and end:
                fragments.append({"start": start, "end": end})
    return fragments
