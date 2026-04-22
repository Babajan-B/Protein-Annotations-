from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from html import escape
import math
import re
from typing import Any
from urllib.parse import quote

import requests

from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import MissenseVariant


PRIMATE_PRIORITY = [
    "homo_sapiens",
    "pan_troglodytes",
    "pan_paniscus",
    "gorilla_gorilla",
    "pongo_abelii",
    "nomascus_leucogenys",
    "macaca_mulatta",
    "papio_anubis",
    "chlorocebus_sabaeus",
    "callithrix_jacchus",
]

DISPLAY_SPECIES = {
    "homo_sapiens": "Human",
    "pan_troglodytes": "Chimpanzee",
    "pan_paniscus": "Bonobo",
    "gorilla_gorilla": "Gorilla",
    "pongo_abelii": "Sumatran orangutan",
    "nomascus_leucogenys": "Northern white-cheeked gibbon",
    "macaca_mulatta": "Rhesus macaque",
    "papio_anubis": "Olive baboon",
    "chlorocebus_sabaeus": "Green monkey",
    "callithrix_jacchus": "Common marmoset",
}

GRCH38_CHROMOSOME_LENGTHS = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "MT": 16569,
    "M": 16569,
}


@dataclass(frozen=True)
class InsightNotice:
    title: str
    body: str


@dataclass
class _TreeNode:
    members: list[str]
    left: "_TreeNode | None" = None
    right: "_TreeNode | None" = None
    height: float = 0.0
    name: str | None = None
    x: float = 0.0
    y: float = 0.0
    angle: float = 0.0


class GenomicInsightService:
    def __init__(
        self,
        *,
        ensembl_base_url: str,
        ncbi_datasets_base_url: str,
        timeout_seconds: float,
        user_agent: str,
    ) -> None:
        self.ensembl_base_url = ensembl_base_url.rstrip("/")
        self.ncbi_datasets_base_url = ncbi_datasets_base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.session = requests.Session()
        self.headers = {"Accept": "application/json", "User-Agent": user_agent}

    def collect(self, protein: ProteinRecord, variant: MissenseVariant) -> dict[str, Any]:
        bundle: dict[str, Any] = {
            "conservation_panel": None,
            "phylogenetic_panel": None,
            "chromosome_location_panel": None,
            "insight_notices": [],
        }

        symbol = protein.gene_symbol or protein.entry_name
        if not symbol:
            bundle["insight_notices"].append(
                InsightNotice("External context unavailable", "Resolved protein does not expose a stable human gene symbol.").__dict__
            )
            return bundle

        try:
            gene_payload = self._request_json(
                f"/lookup/symbol/homo_sapiens/{quote(symbol)}",
                params={"expand": 1},
            )
        except requests.RequestException as exc:
            bundle["insight_notices"].append(
                InsightNotice("Ensembl gene lookup failed", str(exc)).__dict__
            )
            gene_payload = None

        if gene_payload is not None:
            try:
                chromosome_panel = self._build_chromosome_location_panel(gene_payload, protein, variant)
                bundle["chromosome_location_panel"] = chromosome_panel
            except requests.RequestException as exc:
                bundle["insight_notices"].append(
                    InsightNotice("Ensembl chromosome mapping incomplete", str(exc)).__dict__
                )

        if bundle["chromosome_location_panel"] is None:
            try:
                chromosome_panel = self._build_chromosome_location_panel_ncbi(symbol, protein, variant)
                bundle["chromosome_location_panel"] = chromosome_panel
                bundle["insight_notices"].append(
                    InsightNotice(
                        "Chromosome panel source",
                        "Chromosome and exon mapping is using NCBI Datasets fallback for this report.",
                    ).__dict__
                )
            except requests.RequestException as exc:
                bundle["insight_notices"].append(
                    InsightNotice("Chromosome mapping unavailable", str(exc)).__dict__
                )

        try:
            conservation_panel = self._build_conservation_panel(symbol, protein, variant)
            bundle["conservation_panel"] = conservation_panel
            if conservation_panel is not None:
                bundle["phylogenetic_panel"] = self._build_phylogenetic_panel(conservation_panel["species_rows"])
        except requests.RequestException as exc:
            bundle["insight_notices"].append(
                InsightNotice("Primate conservation unavailable", str(exc)).__dict__
            )

        return bundle

    def _build_chromosome_location_panel(
        self, gene_payload: dict[str, Any], protein: ProteinRecord, variant: MissenseVariant
    ) -> dict[str, Any]:
        transcript = self._select_transcript(gene_payload, protein.length)
        if transcript is None:
            raise requests.RequestException("Could not identify a protein-coding transcript aligned to the resolved sequence.")

        translation = transcript.get("Translation") or transcript.get("translation") or {}
        translation_id = translation.get("id") or translation.get("stable_id") or ""
        if not translation_id:
            raise requests.RequestException("Selected transcript does not include a translation identifier.")

        chromosome = str(
            gene_payload.get("seq_region_name")
            or gene_payload.get("seq_region")
            or transcript.get("seq_region_name")
            or transcript.get("seq_region")
            or ""
        )
        gene_start = int(gene_payload.get("start") or transcript.get("start") or 0)
        gene_end = int(gene_payload.get("end") or transcript.get("end") or 0)
        transcript_start = int(transcript.get("start") or gene_start)
        transcript_end = int(transcript.get("end") or gene_end)
        strand = int(transcript.get("strand") or gene_payload.get("strand") or 1)

        assembly_payload = self._request_json(f"/info/assembly/homo_sapiens/{quote(chromosome)}")
        chromosome_length = self._chromosome_length(assembly_payload, chromosome)

        variant_mapping = self._request_json(f"/map/translation/{quote(translation_id)}/{variant.position}..{variant.position}")
        variant_segments = self._extract_variant_segments(variant_mapping, chromosome)
        if not variant_segments:
            raise requests.RequestException("Protein-to-genome mapping did not return genomic coordinates for the submitted residue.")

        variant_start = min(segment["start"] for segment in variant_segments)
        variant_end = max(segment["end"] for segment in variant_segments)

        exons = transcript.get("Exon") or transcript.get("exons") or []
        if not exons:
            transcript_lookup = self._request_json(f"/lookup/id/{quote(str(transcript['id']))}", params={"expand": 1})
            exons = transcript_lookup.get("Exon") or transcript_lookup.get("exons") or []

        exon_rows = self._compressed_exon_rows(exons, variant_start, variant_end)
        variant_exon = next((row for row in exon_rows if row["contains_variant"]), None)
        exon_focus = self._focus_exon_panel(variant_exon, variant_start, variant_end)

        gene_left_pct = round(((gene_start - 1) / chromosome_length) * 100.0, 4)
        gene_width_pct = max(round(((gene_end - gene_start + 1) / chromosome_length) * 100.0, 4), 0.2)
        variant_left_pct = round(((variant_start - 1) / chromosome_length) * 100.0, 4)

        return {
            "gene_symbol": gene_payload.get("display_name") or protein.gene_symbol or protein.entry_name,
            "gene_id": gene_payload.get("id", ""),
            "transcript_id": transcript.get("id", ""),
            "translation_id": translation_id,
            "chromosome": chromosome,
            "strand_label": "forward" if strand >= 0 else "reverse",
            "chromosome_length": chromosome_length,
            "chromosome_length_label": f"{chromosome_length:,} bp",
            "gene_span_label": f"{gene_start:,}-{gene_end:,}",
            "variant_span_label": f"{chromosome}:{variant_start:,}-{variant_end:,}",
            "gene_left_pct": gene_left_pct,
            "gene_width_pct": gene_width_pct,
            "variant_left_pct": variant_left_pct,
            "gene_track_variant_pct": self._variant_track_pct(transcript_start, transcript_end, variant_start),
            "gene_track_variant_width_pct": self._variant_track_width_pct(transcript_start, transcript_end, variant_start, variant_end),
            "transcript_start": transcript_start,
            "transcript_end": transcript_end,
            "transcript_span_label": f"{transcript_start:,}-{transcript_end:,}",
            "exon_rows": exon_rows,
            "exon_focus": exon_focus,
        }

    def _build_conservation_panel(
        self, symbol: str, protein: ProteinRecord, variant: MissenseVariant
    ) -> dict[str, Any]:
        payload = self._request_json(
            f"/homology/symbol/human/{quote(symbol)}",
            params={
                "type": "orthologues",
                "target_taxon": 9443,
                "sequence": "protein",
                "aligned": 1,
            },
        )

        homology_rows = self._extract_homology_rows(payload, protein.sequence, variant)
        if len(homology_rows) < 2:
            raise requests.RequestException("Not enough primate orthologue rows were returned for a conservation summary.")

        conserved_rows = [row for row in homology_rows if row["site_residue"] not in {"-", "?"}]
        same_as_human = sum(1 for row in conserved_rows if row["site_residue"] == variant.wild_type)
        dominant_residue = Counter(row["site_residue"] for row in conserved_rows).most_common(1)[0][0]
        same_as_mutant = sum(1 for row in conserved_rows if row["site_residue"] == variant.mutant)

        return {
            "human_residue": variant.wild_type,
            "mutant_residue": variant.mutant,
            "variant_position": variant.position,
            "species_rows": homology_rows,
            "orthologue_count": len(homology_rows) - 1,
            "conserved_count": same_as_human,
            "comparable_count": len(conserved_rows),
            "dominant_residue": dominant_residue,
            "mutant_match_count": same_as_mutant,
            "conservation_pct": round((same_as_human / len(conserved_rows)) * 100.0, 1) if conserved_rows else 0.0,
        }

    def _build_phylogenetic_panel(self, species_rows: list[dict[str, Any]]) -> dict[str, Any]:
        sequences = {row["species_label"]: row["sequence"] for row in species_rows if row.get("sequence")}
        root = self._upgma_tree(sequences)
        ordered_labels = [row["species_label"] for row in species_rows]
        matrix_rows = []
        for left_label in ordered_labels:
            left_sequence = sequences[left_label]
            cells = []
            for right_label in ordered_labels:
                identity_pct = round(self._pairwise_identity(left_sequence, sequences[right_label]) * 100.0, 1)
                cells.append(
                    {
                        "species_label": right_label,
                        "value": identity_pct,
                        "class_name": self._matrix_class(identity_pct, left_label == right_label),
                        "is_self": left_label == right_label,
                    }
                )
            matrix_rows.append({"species_label": left_label, "cells": cells})

        non_human_rows = [row for row in species_rows if row["species_label"] != DISPLAY_SPECIES["homo_sapiens"]]
        closest_species = max(non_human_rows, key=lambda row: row["percent_identity"]) if non_human_rows else None
        farthest_species = min(non_human_rows, key=lambda row: row["percent_identity"]) if non_human_rows else None

        return {
            "method": "Whole-gene phylogeny approximated from full-length orthologue protein sequence identity using UPGMA clustering.",
            "species_count": len(species_rows),
            "closest_species": closest_species,
            "farthest_species": farthest_species,
            "matrix_headers": ordered_labels,
            "matrix_rows": matrix_rows,
            "tree_svg": self._build_tree_svg(species_rows),
            "tree_views": self._build_tree_views(species_rows),
            "newick": self._newick_from_tree(root) + ";",
        }

    def _build_tree_views(self, species_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
        return [
            {
                "key": "publication",
                "label": "Publication",
                "summary": "Detailed rectangular tree with residue and identity badges.",
                "svg": self._build_tree_svg(species_rows),
            },
            {
                "key": "compact",
                "label": "Compact",
                "summary": "Dense rectangular tree optimized for quick scanning.",
                "svg": self._build_compact_tree_svg(species_rows),
            },
            {
                "key": "diagonal",
                "label": "Diagonal",
                "summary": "Slanted cladogram view for branch-shape comparison.",
                "svg": self._build_diagonal_tree_svg(species_rows),
            },
            {
                "key": "radial",
                "label": "Radial",
                "summary": "Circular whole-gene tree for visualizing overall divergence.",
                "svg": self._build_radial_tree_svg(species_rows),
            },
        ]

    def _build_chromosome_location_panel_ncbi(
        self, symbol: str, protein: ProteinRecord, variant: MissenseVariant
    ) -> dict[str, Any]:
        payload = self._request_ncbi_json(f"/gene/symbol/{quote(symbol)}/taxon/9606/product_report")
        gene_record = self._extract_ncbi_gene_record(payload)
        transcripts = gene_record.get("transcripts") or []
        transcript = self._select_ncbi_transcript(transcripts, protein)
        if transcript is None:
            raise requests.RequestException("NCBI product report did not expose a protein-coding transcript aligned to the resolved protein.")

        protein_info = transcript.get("protein") or {}
        transcript_id = transcript.get("ensemblTranscript") or transcript.get("accessionVersion") or ""
        translation_id = protein_info.get("ensemblProtein") or protein_info.get("accessionVersion") or ""
        genomic_location = self._select_ncbi_genomic_location(transcript)
        chromosome = self._ncbi_chromosome_label(genomic_location)
        chromosome_length = GRCH38_CHROMOSOME_LENGTHS.get(chromosome)
        if chromosome_length is None:
            raise requests.RequestException(f"Unsupported chromosome label from NCBI report: {chromosome}")

        genomic_range = genomic_location.get("genomicRange") or {}
        gene_start = int(gene_record.get("genomicRanges", [{}])[0].get("range", [{}])[0].get("begin") or genomic_range.get("begin") or 0)
        gene_end = int(gene_record.get("genomicRanges", [{}])[0].get("range", [{}])[0].get("end") or genomic_range.get("end") or 0)
        transcript_start = int(genomic_range.get("begin") or gene_start)
        transcript_end = int(genomic_range.get("end") or gene_end)
        orientation = str(genomic_range.get("orientation") or "plus").lower()
        strand_label = "forward" if orientation == "plus" else "reverse"

        cds_ranges = (((transcript.get("cds") or {}).get("range")) or [])
        cds_start_tx, cds_end_tx = self._select_ncbi_cds_range(cds_ranges)
        codon_start_tx = cds_start_tx + ((variant.position - 1) * 3)
        codon_end_tx = codon_start_tx + 2
        if codon_end_tx > cds_end_tx:
            raise requests.RequestException("NCBI transcript CDS range does not include the submitted protein residue.")

        exons = genomic_location.get("exons") or []
        variant_segments = self._map_ncbi_transcript_interval_to_genome(exons, codon_start_tx, codon_end_tx)
        if not variant_segments:
            raise requests.RequestException("NCBI exon mapping could not localize the submitted residue on the genome.")

        variant_start = min(segment["start"] for segment in variant_segments)
        variant_end = max(segment["end"] for segment in variant_segments)
        exon_rows = self._compressed_exon_rows(exons, variant_start, variant_end)
        variant_exon = next((row for row in exon_rows if row["contains_variant"]), None)
        exon_focus = self._focus_exon_panel(variant_exon, variant_start, variant_end)

        gene_left_pct = round(((gene_start - 1) / chromosome_length) * 100.0, 4)
        gene_width_pct = max(round(((gene_end - gene_start + 1) / chromosome_length) * 100.0, 4), 0.2)
        variant_left_pct = round(((variant_start - 1) / chromosome_length) * 100.0, 4)

        return {
            "gene_symbol": gene_record.get("symbol") or protein.gene_symbol or protein.entry_name,
            "gene_id": str(gene_record.get("geneId") or ""),
            "transcript_id": transcript_id,
            "translation_id": translation_id,
            "chromosome": chromosome,
            "strand_label": strand_label,
            "chromosome_length": chromosome_length,
            "chromosome_length_label": f"{chromosome_length:,} bp",
            "gene_span_label": f"{gene_start:,}-{gene_end:,}",
            "variant_span_label": f"{chromosome}:{variant_start:,}-{variant_end:,}",
            "gene_left_pct": gene_left_pct,
            "gene_width_pct": gene_width_pct,
            "variant_left_pct": variant_left_pct,
            "gene_track_variant_pct": self._variant_track_pct(transcript_start, transcript_end, variant_start),
            "gene_track_variant_width_pct": self._variant_track_width_pct(transcript_start, transcript_end, variant_start, variant_end),
            "transcript_start": transcript_start,
            "transcript_end": transcript_end,
            "transcript_span_label": f"{transcript_start:,}-{transcript_end:,}",
            "exon_rows": exon_rows,
            "exon_focus": exon_focus,
            "source": "NCBI Datasets",
        }

    def _extract_homology_rows(
        self, payload: dict[str, Any], human_sequence: str, variant: MissenseVariant
    ) -> list[dict[str, Any]]:
        data_rows = payload.get("data") or []
        homologies = data_rows[0].get("homologies", []) if data_rows else []
        rows = [
            {
                "species_key": "homo_sapiens",
                "species_label": DISPLAY_SPECIES["homo_sapiens"],
                "site_residue": variant.wild_type,
                "sequence": human_sequence,
                "percent_identity": 100.0,
                "match_class": "match",
            }
        ]
        seen = {"homo_sapiens"}

        for priority_key in PRIMATE_PRIORITY[1:]:
            for homology in homologies:
                target = homology.get("target", {})
                species_key = str(target.get("species") or "")
                if species_key != priority_key or species_key in seen:
                    continue

                source_align = str((homology.get("source") or {}).get("align_seq") or "")
                target_align = str(target.get("align_seq") or "")
                target_sequence = target_align.replace("-", "") or str(target.get("seq") or "")
                site_residue = self._mapped_alignment_residue(source_align, target_align, variant.position)
                percent_identity = float(target.get("perc_id") or homology.get("target_perc_id") or 0.0)

                rows.append(
                    {
                        "species_key": species_key,
                        "species_label": DISPLAY_SPECIES.get(species_key, species_key.replace("_", " ").title()),
                        "site_residue": site_residue or "?",
                        "sequence": target_sequence,
                        "percent_identity": round(percent_identity, 2),
                        "match_class": self._match_class(site_residue, variant.wild_type),
                    }
                )
                seen.add(species_key)
                break

        return rows

    @staticmethod
    def _mapped_alignment_residue(source_align: str, target_align: str, protein_position: int) -> str:
        residue_index = 0
        for source_residue, target_residue in zip(source_align, target_align):
            if source_residue != "-":
                residue_index += 1
            if residue_index == protein_position:
                return target_residue
        return "?"

    @staticmethod
    def _match_class(site_residue: str, reference_residue: str) -> str:
        if site_residue == reference_residue:
            return "match"
        if site_residue in {"-", "?"}:
            return "unknown"
        return "mismatch"

    @staticmethod
    def _extract_ncbi_gene_record(payload: dict[str, Any]) -> dict[str, Any]:
        for key in ("reports", "genes", "geneDescriptors", "gene_descriptors"):
            value = payload.get(key)
            if isinstance(value, list) and value:
                return value[0]
        if payload.get("geneId") or payload.get("symbol"):
            return payload
        raise requests.RequestException("NCBI product report did not contain a gene record.")

    @staticmethod
    def _select_ncbi_transcript(transcripts: list[dict[str, Any]], protein: ProteinRecord) -> dict[str, Any] | None:
        refseq_targets = {accession.split(".")[0] for accession in protein.refseq_proteins}
        candidates: list[tuple[tuple[int, int, int, int], dict[str, Any]]] = []
        for transcript in transcripts:
            protein_info = transcript.get("protein") or {}
            protein_accession = str(protein_info.get("accessionVersion") or "").split(".")[0]
            protein_length = int(protein_info.get("length") or 0)
            accession_flag = 0 if protein_accession and protein_accession in refseq_targets else 1
            length_flag = 0 if protein_length == protein.length else 1
            type_flag = 0 if str(transcript.get("type") or "").upper() == "PROTEIN_CODING" else 1
            score = (accession_flag, length_flag, type_flag, abs(protein_length - protein.length))
            candidates.append((score, transcript))

        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0])
        return candidates[0][1]

    @staticmethod
    def _select_ncbi_genomic_location(transcript: dict[str, Any]) -> dict[str, Any]:
        locations = transcript.get("genomicLocations") or []
        if not locations:
            raise requests.RequestException("NCBI transcript record did not include genomic locations.")
        primary = next(
            (
                location
                for location in locations
                if "Primary Assembly" in str(location.get("sequenceName") or "")
            ),
            None,
        )
        return primary or locations[0]

    @staticmethod
    def _select_ncbi_cds_range(cds_ranges: list[dict[str, Any]]) -> tuple[int, int]:
        if not cds_ranges:
            raise requests.RequestException("NCBI transcript record did not include CDS ranges.")
        begins = [int(item["begin"]) for item in cds_ranges if item.get("begin") is not None]
        ends = [int(item["end"]) for item in cds_ranges if item.get("end") is not None]
        if not begins or not ends:
            raise requests.RequestException("NCBI transcript CDS range was incomplete.")
        return min(begins), max(ends)

    @staticmethod
    def _ncbi_chromosome_label(genomic_location: dict[str, Any]) -> str:
        sequence_name = str(genomic_location.get("sequenceName") or "")
        match = re.search(r"Chromosome\s+([0-9XYM]+)", sequence_name, re.IGNORECASE)
        if match:
            return match.group(1).upper()
        accession = str(genomic_location.get("genomicAccessionVersion") or "")
        accession_match = re.match(r"NC_0+([0-9]+)\.", accession)
        if accession_match:
            number = int(accession_match.group(1))
            if number == 23:
                return "X"
            if number == 24:
                return "Y"
            return str(number)
        raise requests.RequestException("Could not derive a chromosome label from the NCBI genomic location.")

    @staticmethod
    def _map_ncbi_transcript_interval_to_genome(
        exons: list[dict[str, Any]], transcript_start: int, transcript_end: int
    ) -> list[dict[str, int]]:
        ordered_exons = sorted(exons, key=lambda exon: int(exon.get("order") or 0))
        remaining_start = transcript_start
        remaining_end = transcript_end
        transcript_cursor = 1
        segments: list[dict[str, int]] = []

        for exon in ordered_exons:
            exon_begin = int(exon["begin"])
            exon_end = int(exon["end"])
            exon_length = abs(exon_end - exon_begin) + 1
            exon_tx_start = transcript_cursor
            exon_tx_end = transcript_cursor + exon_length - 1

            overlap_start = max(remaining_start, exon_tx_start)
            overlap_end = min(remaining_end, exon_tx_end)
            if overlap_start <= overlap_end:
                offset_start = overlap_start - exon_tx_start
                offset_end = overlap_end - exon_tx_start
                orientation = str(exon.get("orientation") or "plus").lower()
                if orientation == "minus":
                    genomic_segment_start = exon_end - offset_end
                    genomic_segment_end = exon_end - offset_start
                else:
                    genomic_segment_start = exon_begin + offset_start
                    genomic_segment_end = exon_begin + offset_end
                segments.append(
                    {
                        "start": min(genomic_segment_start, genomic_segment_end),
                        "end": max(genomic_segment_start, genomic_segment_end),
                    }
                )

            transcript_cursor = exon_tx_end + 1
            if transcript_cursor > remaining_end:
                break

        return segments

    @staticmethod
    def _chromosome_length(payload: dict[str, Any], chromosome: str) -> int:
        if payload.get("length"):
            return int(payload["length"])
        for region in payload.get("top_level_region", []):
            if str(region.get("name")) == chromosome:
                return int(region["length"])
        raise requests.RequestException(f"Chromosome length not found for {chromosome}.")

    @staticmethod
    def _extract_variant_segments(payload: dict[str, Any], chromosome: str) -> list[dict[str, int]]:
        segments = []
        for mapping in payload.get("mappings", []):
            mapped = mapping.get("mapped") or mapping
            seq_region = str(mapped.get("seq_region_name") or mapped.get("seq_region") or "")
            if chromosome and seq_region and seq_region != chromosome:
                continue
            start = mapped.get("start")
            end = mapped.get("end")
            if start is None or end is None:
                continue
            segments.append({"start": int(start), "end": int(end)})
        return segments

    @staticmethod
    def _select_transcript(gene_payload: dict[str, Any], protein_length: int) -> dict[str, Any] | None:
        transcripts = gene_payload.get("Transcript") or gene_payload.get("transcripts") or []
        canonical_transcript = gene_payload.get("canonical_transcript") or gene_payload.get("canonicalTranscript") or ""

        candidates: list[tuple[tuple[int, int, int, int], dict[str, Any]]] = []
        for transcript in transcripts:
            translation = transcript.get("Translation") or transcript.get("translation")
            if not translation:
                continue
            transcript_id = str(transcript.get("id") or transcript.get("stable_id") or "")
            translation_length = int(translation.get("length") or translation.get("Length") or 0)
            canonical_flag = 0 if transcript.get("is_canonical") or transcript_id == canonical_transcript else 1
            length_flag = 0 if translation_length == protein_length else 1
            score = (length_flag, abs(translation_length - protein_length), canonical_flag, -translation_length)
            candidates.append((score, transcript))

        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0])
        return candidates[0][1]

    @staticmethod
    def _compressed_exon_rows(
        exons: list[dict[str, Any]], variant_start: int, variant_end: int
    ) -> list[dict[str, Any]]:
        normalized = []
        sortable_exons = sorted(
            exons,
            key=lambda row: int(row.get("start") or row.get("begin") or 0),
        )
        for index, exon in enumerate(sortable_exons, start=1):
            start = int(exon.get("start") or exon.get("begin") or 0)
            end = int(exon.get("end") or exon.get("stop") or exon.get("begin") or 0)
            normalized.append(
                {
                    "number": int(exon.get("order") or index),
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "contains_variant": not (variant_end < start or variant_start > end),
                }
            )

        display_cursor = 0.0
        for index, exon in enumerate(normalized):
            exon_display = max(float(exon["length"]), 180.0)
            exon["display_start"] = display_cursor
            exon["display_width"] = exon_display
            display_cursor += exon_display
            if index < len(normalized) - 1:
                intron = normalized[index + 1]["start"] - exon["end"] - 1
                display_cursor += min(max(float(intron), 24.0), 160.0)

        total_display = display_cursor or 1.0
        rows = []
        for exon in normalized:
            row = {
                **exon,
                "left_pct": round((float(exon["display_start"]) / total_display) * 100.0, 3),
                "width_pct": round((float(exon["display_width"]) / total_display) * 100.0, 3),
            }
            rows.append(row)
        return rows

    @staticmethod
    def _focus_exon_panel(
        exon_row: dict[str, Any] | None, variant_start: int, variant_end: int
    ) -> dict[str, Any] | None:
        if exon_row is None:
            return None

        exon_start = int(exon_row["start"])
        exon_end = int(exon_row["end"])
        exon_length = exon_end - exon_start + 1
        flank = 120
        visible_start = exon_start
        visible_end = exon_end
        left_truncated = False
        right_truncated = False

        if exon_length > (flank * 2) + 30:
            visible_start = max(exon_start, variant_start - flank)
            visible_end = min(exon_end, variant_end + flank)
            left_truncated = visible_start > exon_start
            right_truncated = visible_end < exon_end

        visible_length = max(visible_end - visible_start + 1, 1)
        variant_left_pct = round(((variant_start - visible_start) / visible_length) * 100.0, 3)
        variant_width_pct = max(round(((variant_end - variant_start + 1) / visible_length) * 100.0, 3), 1.1)

        return {
            "label": f"Exon {exon_row['number']}",
            "start": exon_start,
            "end": exon_end,
            "length": exon_length,
            "visible_start": visible_start,
            "visible_end": visible_end,
            "left_truncated": left_truncated,
            "right_truncated": right_truncated,
            "variant_left_pct": variant_left_pct,
            "variant_width_pct": variant_width_pct,
        }

    @staticmethod
    def _variant_track_pct(track_start: int, track_end: int, variant_start: int) -> float:
        span = max(track_end - track_start + 1, 1)
        return round(((variant_start - track_start) / span) * 100.0, 3)

    @staticmethod
    def _variant_track_width_pct(track_start: int, track_end: int, variant_start: int, variant_end: int) -> float:
        span = max(track_end - track_start + 1, 1)
        return max(round(((variant_end - variant_start + 1) / span) * 100.0, 3), 0.8)

    def _request_json(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        url = f"{self.ensembl_base_url}{path}"
        try:
            response = self.session.get(url, params=params, headers=self.headers, timeout=self.timeout_seconds)
            response.raise_for_status()
        except requests.RequestException as exc:
            raise requests.RequestException(f"Ensembl request failed for {path}: {exc}") from exc
        return response.json()

    def _request_ncbi_json(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        url = f"{self.ncbi_datasets_base_url}{path}"
        try:
            response = self.session.get(url, params=params, headers=self.headers, timeout=self.timeout_seconds)
            response.raise_for_status()
        except requests.RequestException as exc:
            raise requests.RequestException(f"NCBI Datasets request failed for {path}: {exc}") from exc
        return response.json()

    @staticmethod
    def _tree_metadata(species_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        return {
            row["species_label"]: {
                "site_residue": row.get("site_residue", "?"),
                "percent_identity": float(row.get("percent_identity", 0.0)),
                "match_class": row.get("match_class", "unknown"),
            }
            for row in species_rows
        }

    @staticmethod
    def _tree_sequences(species_rows: list[dict[str, Any]]) -> dict[str, str]:
        return {row["species_label"]: row["sequence"] for row in species_rows if row.get("sequence")}

    def _build_tree_svg(self, species_rows: list[dict[str, Any]]) -> str:
        metadata = self._tree_metadata(species_rows)
        sequences = self._tree_sequences(species_rows)
        if len(sequences) == 1:
            label = next(iter(sequences))
            return (
                '<svg viewBox="0 0 520 80" class="phylo-svg" role="img" aria-label="Single-species tree">'
                f'<text x="24" y="42" font-size="14" fill="#1e2b26">{escape(label)}</text>'
                "</svg>"
            )

        root = self._upgma_tree(sequences)
        self._assign_tree_coordinates(root)

        leaf_count = len(sequences)
        height = max((leaf_count * 50) + 110, 230)
        width = 980
        branch_width = 360
        max_height = max(root.height, 0.01)
        self._scale_tree_x(root, max_height=max_height, branch_width=branch_width)

        segments: list[str] = []
        labels: list[str] = []
        guides: list[str] = []
        leaf_positions = self._leaf_positions(root)
        label_column_x = branch_width + 72
        residue_column_x = width - 210
        identity_column_x = width - 108

        for species_label in self._leaf_order(root):
            y = leaf_positions[species_label]
            highlight_fill = "#f0f7f3" if species_label == "Human" else "transparent"
            guides.append(
                f'<rect x="16" y="{y - 18:.1f}" width="{width - 32}" height="36" rx="12" fill="{highlight_fill}"></rect>'
            )
            guides.append(
                f'<line x1="{label_column_x - 16:.1f}" y1="{y:.1f}" x2="{residue_column_x - 20:.1f}" y2="{y:.1f}" '
                'stroke="#d6e1da" stroke-width="1" stroke-dasharray="4 6"></line>'
            )

        self._render_tree(
            root,
            segments,
            labels,
            metadata=metadata,
            label_column_x=label_column_x,
            residue_column_x=residue_column_x,
            identity_column_x=identity_column_x,
        )

        header = (
            '<text x="24" y="28" font-size="16" font-weight="700" fill="#15352d">Whole-gene primate phylogeny</text>'
            '<text x="24" y="48" font-size="11" fill="#5d7268">Leaf labels carry the residue at the submitted site and the full-length identity to human.</text>'
            f'<text x="{label_column_x:.1f}" y="76" font-size="11" font-weight="700" fill="#5d7268">Species</text>'
            f'<text x="{residue_column_x:.1f}" y="76" font-size="11" font-weight="700" fill="#5d7268">Site residue</text>'
            f'<text x="{identity_column_x:.1f}" y="76" font-size="11" font-weight="700" fill="#5d7268">Gene identity</text>'
        )
        scale_bar = self._tree_scale_bar(max_height=max_height, x=36, y=height - 22, width=120)

        return (
            f'<svg viewBox="0 0 {width} {height}" class="phylo-svg" role="img" aria-label="Primate protein tree">'
            '<defs>'
            '<linearGradient id="phylo-bg" x1="0%" y1="0%" x2="100%" y2="100%">'
            '<stop offset="0%" stop-color="#fbfcf9"></stop>'
            '<stop offset="100%" stop-color="#f1f6f3"></stop>'
            '</linearGradient>'
            '</defs>'
            '<rect x="0" y="0" width="100%" height="100%" rx="24" fill="url(#phylo-bg)"></rect>'
            f"{header}"
            f'{"".join(guides)}'
            f'{"".join(segments)}'
            f'{"".join(labels)}'
            f"{scale_bar}"
            "</svg>"
        )

    def _build_compact_tree_svg(self, species_rows: list[dict[str, Any]]) -> str:
        metadata = self._tree_metadata(species_rows)
        sequences = self._tree_sequences(species_rows)
        if len(sequences) == 1:
            label = next(iter(sequences))
            return (
                '<svg viewBox="0 0 420 80" class="phylo-svg" role="img" aria-label="Single-species compact tree">'
                f'<text x="24" y="42" font-size="14" fill="#1e2b26">{escape(label)}</text>'
                "</svg>"
            )

        root = self._upgma_tree(sequences)
        self._assign_tree_coordinates(root, start_y=42.0, spacing=30.0)
        leaf_count = len(sequences)
        height = max((leaf_count * 34) + 82, 180)
        width = 720
        branch_width = 230
        max_height = max(root.height, 0.01)
        self._scale_tree_x(root, max_height=max_height, branch_width=branch_width, start_x=30.0)

        segments: list[str] = []
        labels: list[str] = []
        self._render_compact_tree(root, segments, labels, metadata=metadata, label_x=branch_width + 44)
        return (
            f'<svg viewBox="0 0 {width} {height}" class="phylo-svg" role="img" aria-label="Compact primate tree">'
            '<rect x="0" y="0" width="100%" height="100%" rx="20" fill="#f8faf7"></rect>'
            '<text x="24" y="24" font-size="14" font-weight="700" fill="#15352d">Compact whole-gene tree</text>'
            f'{"".join(segments)}'
            f'{"".join(labels)}'
            "</svg>"
        )

    def _build_diagonal_tree_svg(self, species_rows: list[dict[str, Any]]) -> str:
        metadata = self._tree_metadata(species_rows)
        sequences = self._tree_sequences(species_rows)
        if len(sequences) == 1:
            label = next(iter(sequences))
            return (
                '<svg viewBox="0 0 420 80" class="phylo-svg" role="img" aria-label="Single-species diagonal tree">'
                f'<text x="24" y="42" font-size="14" fill="#1e2b26">{escape(label)}</text>'
                "</svg>"
            )

        root = self._upgma_tree(sequences)
        self._assign_tree_coordinates(root, start_y=56.0, spacing=34.0)
        leaf_count = len(sequences)
        height = max((leaf_count * 40) + 92, 220)
        width = 760
        branch_width = 270
        max_height = max(root.height, 0.01)
        self._scale_tree_x(root, max_height=max_height, branch_width=branch_width, start_x=44.0)

        segments: list[str] = []
        labels: list[str] = []
        self._render_diagonal_tree(root, segments, labels, metadata=metadata, label_x=branch_width + 52)
        return (
            f'<svg viewBox="0 0 {width} {height}" class="phylo-svg" role="img" aria-label="Diagonal primate tree">'
            '<defs><linearGradient id="diag-bg" x1="0%" y1="0%" x2="100%" y2="100%"><stop offset="0%" stop-color="#fbfcf9"></stop><stop offset="100%" stop-color="#eef5f1"></stop></linearGradient></defs>'
            '<rect x="0" y="0" width="100%" height="100%" rx="20" fill="url(#diag-bg)"></rect>'
            '<text x="24" y="28" font-size="14" font-weight="700" fill="#15352d">Diagonal cladogram</text>'
            f'{"".join(segments)}'
            f'{"".join(labels)}'
            "</svg>"
        )

    def _build_radial_tree_svg(self, species_rows: list[dict[str, Any]]) -> str:
        metadata = self._tree_metadata(species_rows)
        sequences = self._tree_sequences(species_rows)
        if len(sequences) == 1:
            label = next(iter(sequences))
            return (
                '<svg viewBox="0 0 420 120" class="phylo-svg" role="img" aria-label="Single-species radial tree">'
                f'<circle cx="80" cy="60" r="8" fill="#114f43"></circle><text x="100" y="65" font-size="14" fill="#1e2b26">{escape(label)}</text>'
                "</svg>"
            )

        root = self._upgma_tree(sequences)
        center_x = 340.0
        center_y = 300.0
        inner_radius = 44.0
        outer_radius = 220.0
        max_height = max(root.height, 0.01)
        self._assign_radial_angles(root, start_angle=-1.45, end_angle=1.45)
        self._set_radial_layout(root, center_x=center_x, center_y=center_y, inner_radius=inner_radius, outer_radius=outer_radius, max_height=max_height)

        segments: list[str] = []
        labels: list[str] = []
        self._render_radial_tree(root, segments, labels, metadata=metadata, center_x=center_x, center_y=center_y)
        legend_items = [
            '<rect x="580" y="90" width="22" height="22" rx="11" fill="#d6efe3"></rect><text x="612" y="105" font-size="12" fill="#1e2b26">Site matches human</text>',
            '<rect x="580" y="122" width="22" height="22" rx="11" fill="#f5d8d0"></rect><text x="612" y="137" font-size="12" fill="#1e2b26">Site differs from human</text>',
            '<rect x="580" y="154" width="22" height="22" rx="11" fill="#114f43"></rect><text x="612" y="169" font-size="12" fill="#1e2b26">Human reference lineage</text>',
        ]
        return (
            '<svg viewBox="0 0 920 620" class="phylo-svg" role="img" aria-label="Radial primate tree">'
            '<defs><radialGradient id="radial-bg" cx="45%" cy="45%" r="70%"><stop offset="0%" stop-color="#fbfcf9"></stop><stop offset="100%" stop-color="#eef5f1"></stop></radialGradient></defs>'
            '<rect x="0" y="0" width="100%" height="100%" rx="24" fill="url(#radial-bg)"></rect>'
            '<text x="26" y="32" font-size="14" font-weight="700" fill="#15352d">Radial whole-gene tree</text>'
            f'{"".join(segments)}'
            f'{"".join(labels)}'
            f'{"".join(legend_items)}'
            "</svg>"
        )

    def _upgma_tree(self, sequences: dict[str, str]) -> _TreeNode:
        clusters = {
            name: _TreeNode(members=[name], name=name, height=0.0)
            for name in sequences
        }
        sizes = {name: 1 for name in sequences}
        distances: dict[frozenset[str], float] = {}
        names = list(sequences)
        for index, left in enumerate(names):
            for right in names[index + 1 :]:
                distances[frozenset({left, right})] = self._sequence_distance(sequences[left], sequences[right])

        counter = 0
        while len(clusters) > 1:
            pair = min(distances, key=distances.get)
            left_name, right_name = tuple(pair)
            left_node = clusters.pop(left_name)
            right_node = clusters.pop(right_name)
            left_size = sizes.pop(left_name)
            right_size = sizes.pop(right_name)
            new_name = f"cluster_{counter}"
            counter += 1
            pair_distance = distances[pair]
            new_node = _TreeNode(
                members=left_node.members + right_node.members,
                left=left_node,
                right=right_node,
                height=pair_distance / 2.0,
            )

            keys_to_remove = [key for key in distances if left_name in key or right_name in key]
            for key in keys_to_remove:
                distances.pop(key, None)

            for other_name, other_node in clusters.items():
                other_distance_left = self._cluster_distance(left_node.members, other_node.members, sequences)
                other_distance_right = self._cluster_distance(right_node.members, other_node.members, sequences)
                avg_distance = ((other_distance_left * left_size) + (other_distance_right * right_size)) / (left_size + right_size)
                distances[frozenset({new_name, other_name})] = avg_distance

            clusters[new_name] = new_node
            sizes[new_name] = left_size + right_size

        return next(iter(clusters.values()))

    def _cluster_distance(self, members_a: list[str], members_b: list[str], sequences: dict[str, str]) -> float:
        values = []
        for left in members_a:
            for right in members_b:
                values.append(self._sequence_distance(sequences[left], sequences[right]))
        return sum(values) / len(values)

    @staticmethod
    def _sequence_distance(sequence_a: str, sequence_b: str) -> float:
        min_length = min(len(sequence_a), len(sequence_b))
        if min_length == 0:
            return 1.0
        matches = sum(1 for left, right in zip(sequence_a[:min_length], sequence_b[:min_length]) if left == right)
        identity = matches / min_length
        length_penalty = abs(len(sequence_a) - len(sequence_b)) / max(len(sequence_a), len(sequence_b))
        return max(0.0, 1.0 - identity + (length_penalty * 0.25))

    @staticmethod
    def _pairwise_identity(sequence_a: str, sequence_b: str) -> float:
        min_length = min(len(sequence_a), len(sequence_b))
        if min_length == 0:
            return 0.0
        matches = sum(1 for left, right in zip(sequence_a[:min_length], sequence_b[:min_length]) if left == right)
        return matches / min_length

    @staticmethod
    def _matrix_class(identity_pct: float, is_self: bool) -> str:
        if is_self:
            return "matrix-self"
        if identity_pct >= 97.0:
            return "matrix-vhigh"
        if identity_pct >= 92.0:
            return "matrix-high"
        if identity_pct >= 85.0:
            return "matrix-medium"
        return "matrix-low"

    def _assign_tree_coordinates(self, root: _TreeNode, *, start_y: float = 40.0, spacing: float = 38.0) -> None:
        leaves = self._leaf_order(root)
        y_positions = {name: start_y + (index * spacing) for index, name in enumerate(leaves)}
        self._set_tree_y(root, y_positions)

    def _leaf_order(self, node: _TreeNode) -> list[str]:
        if node.name is not None:
            return [node.name]
        return self._leaf_order(node.left) + self._leaf_order(node.right)

    def _set_tree_y(self, node: _TreeNode, positions: dict[str, float]) -> float:
        if node.name is not None:
            node.y = positions[node.name]
            return node.y
        node.y = (self._set_tree_y(node.left, positions) + self._set_tree_y(node.right, positions)) / 2.0
        return node.y

    def _scale_tree_x(self, node: _TreeNode, *, max_height: float, branch_width: float, start_x: float = 32.0) -> None:
        node.x = start_x + (branch_width * (1.0 - (node.height / max_height)))
        if node.left is not None:
            self._scale_tree_x(node.left, max_height=max_height, branch_width=branch_width, start_x=start_x)
        if node.right is not None:
            self._scale_tree_x(node.right, max_height=max_height, branch_width=branch_width, start_x=start_x)

    def _render_tree(
        self,
        node: _TreeNode,
        segments: list[str],
        labels: list[str],
        *,
        metadata: dict[str, dict[str, Any]],
        label_column_x: float,
        residue_column_x: float,
        identity_column_x: float,
    ) -> None:
        if node.left is not None and node.right is not None:
            segments.append(
                f'<line x1="{node.x:.1f}" y1="{node.left.y:.1f}" x2="{node.x:.1f}" y2="{node.right.y:.1f}" '
                'stroke="#5d7268" stroke-width="2.3"></line>'
            )
            segments.append(
                f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="3.4" fill="#f3efe6" stroke="#5d7268" stroke-width="1.3"></circle>'
            )
            for child in (node.left, node.right):
                segments.append(
                    f'<line x1="{node.x:.1f}" y1="{child.y:.1f}" x2="{child.x:.1f}" y2="{child.y:.1f}" '
                    'stroke="#5d7268" stroke-width="2.3"></line>'
                )
                self._render_tree(
                    child,
                    segments,
                    labels,
                    metadata=metadata,
                    label_column_x=label_column_x,
                    residue_column_x=residue_column_x,
                    identity_column_x=identity_column_x,
                )
            return

        label = escape(node.name or "")
        species_label = node.name or ""
        meta = metadata.get(species_label, {})
        percent_identity = float(meta.get("percent_identity", 0.0))
        identity_fill, identity_text = self._identity_badge_colors(percent_identity, species_label == "Human")
        residue_fill, residue_text = self._residue_badge_colors(meta.get("match_class", "unknown"), species_label == "Human")
        site_residue = escape(str(meta.get("site_residue", "?")))
        label_weight = "700" if species_label == "Human" else "500"
        label_fill = "#114f43" if species_label == "Human" else "#1e2b26"

        labels.append(
            f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="5.3" fill="{residue_fill}" stroke="#ffffff" stroke-width="1.6"></circle>'
        )
        labels.append(
            f'<text x="{label_column_x:.1f}" y="{node.y + 4:.1f}" font-size="13" font-weight="{label_weight}" fill="{label_fill}">{label}</text>'
        )
        labels.append(
            f'<rect x="{residue_column_x:.1f}" y="{node.y - 12:.1f}" width="44" height="24" rx="12" fill="{residue_fill}"></rect>'
            f'<text x="{residue_column_x + 22:.1f}" y="{node.y + 4:.1f}" text-anchor="middle" font-size="12" font-weight="700" fill="{residue_text}">{site_residue}</text>'
        )
        labels.append(
            f'<rect x="{identity_column_x:.1f}" y="{node.y - 12:.1f}" width="68" height="24" rx="12" fill="{identity_fill}"></rect>'
            f'<text x="{identity_column_x + 34:.1f}" y="{node.y + 4:.1f}" text-anchor="middle" font-size="11" font-weight="700" fill="{identity_text}">{percent_identity:.1f}%</text>'
        )

    def _render_compact_tree(
        self,
        node: _TreeNode,
        segments: list[str],
        labels: list[str],
        *,
        metadata: dict[str, dict[str, Any]],
        label_x: float,
    ) -> None:
        if node.left is not None and node.right is not None:
            segments.append(
                f'<line x1="{node.x:.1f}" y1="{node.left.y:.1f}" x2="{node.x:.1f}" y2="{node.right.y:.1f}" stroke="#536960" stroke-width="1.8"></line>'
            )
            for child in (node.left, node.right):
                segments.append(
                    f'<line x1="{node.x:.1f}" y1="{child.y:.1f}" x2="{child.x:.1f}" y2="{child.y:.1f}" stroke="#536960" stroke-width="1.8"></line>'
                )
                self._render_compact_tree(child, segments, labels, metadata=metadata, label_x=label_x)
            return

        species_label = node.name or ""
        meta = metadata.get(species_label, {})
        label = escape(species_label)
        residue_fill, residue_text = self._residue_badge_colors(meta.get("match_class", "unknown"), species_label == "Human")
        label_weight = "700" if species_label == "Human" else "500"
        label_fill = "#114f43" if species_label == "Human" else "#1e2b26"
        labels.append(
            f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="4.8" fill="{residue_fill}" stroke="#ffffff" stroke-width="1.2"></circle>'
            f'<text x="{label_x:.1f}" y="{node.y + 4:.1f}" font-size="12" font-weight="{label_weight}" fill="{label_fill}">{label}</text>'
            f'<text x="{label_x + 166:.1f}" y="{node.y + 4:.1f}" font-size="11" fill="{residue_text}">{escape(str(meta.get("site_residue", "?")))}</text>'
        )

    def _render_diagonal_tree(
        self,
        node: _TreeNode,
        segments: list[str],
        labels: list[str],
        *,
        metadata: dict[str, dict[str, Any]],
        label_x: float,
    ) -> None:
        if node.left is not None and node.right is not None:
            segments.append(
                f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="3.1" fill="#f3efe6" stroke="#5d7268" stroke-width="1.1"></circle>'
            )
            for child in (node.left, node.right):
                segments.append(
                    f'<line x1="{node.x:.1f}" y1="{node.y:.1f}" x2="{child.x:.1f}" y2="{child.y:.1f}" stroke="#5d7268" stroke-width="2"></line>'
                )
                self._render_diagonal_tree(child, segments, labels, metadata=metadata, label_x=label_x)
            return

        species_label = node.name or ""
        meta = metadata.get(species_label, {})
        residue_fill, residue_text = self._residue_badge_colors(meta.get("match_class", "unknown"), species_label == "Human")
        label_weight = "700" if species_label == "Human" else "500"
        label_fill = "#114f43" if species_label == "Human" else "#1e2b26"
        labels.append(
            f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="5.2" fill="{residue_fill}" stroke="#ffffff" stroke-width="1.4"></circle>'
            f'<text x="{label_x:.1f}" y="{node.y + 4:.1f}" font-size="12.5" font-weight="{label_weight}" fill="{label_fill}">{escape(species_label)}</text>'
            f'<text x="{label_x + 164:.1f}" y="{node.y + 4:.1f}" font-size="11" fill="{residue_text}">{float(meta.get("percent_identity", 0.0)):.1f}%</text>'
        )

    def _assign_radial_angles(self, root: _TreeNode, *, start_angle: float, end_angle: float) -> None:
        leaves = self._leaf_order(root)
        if len(leaves) == 1:
            angle_map = {leaves[0]: 0.0}
        else:
            step = (end_angle - start_angle) / max(len(leaves) - 1, 1)
            angle_map = {name: start_angle + (index * step) for index, name in enumerate(leaves)}
        self._set_tree_angle(root, angle_map)

    def _set_tree_angle(self, node: _TreeNode, angle_map: dict[str, float]) -> float:
        if node.name is not None:
            node.angle = angle_map[node.name]
            return node.angle
        node.angle = (self._set_tree_angle(node.left, angle_map) + self._set_tree_angle(node.right, angle_map)) / 2.0
        return node.angle

    def _set_radial_layout(
        self,
        node: _TreeNode,
        *,
        center_x: float,
        center_y: float,
        inner_radius: float,
        outer_radius: float,
        max_height: float,
    ) -> None:
        radius = inner_radius + ((1.0 - (node.height / max_height)) * outer_radius)
        node.x = center_x + (radius * math.cos(node.angle))
        node.y = center_y + (radius * math.sin(node.angle))
        if node.left is not None:
            self._set_radial_layout(
                node.left,
                center_x=center_x,
                center_y=center_y,
                inner_radius=inner_radius,
                outer_radius=outer_radius,
                max_height=max_height,
            )
        if node.right is not None:
            self._set_radial_layout(
                node.right,
                center_x=center_x,
                center_y=center_y,
                inner_radius=inner_radius,
                outer_radius=outer_radius,
                max_height=max_height,
            )

    def _render_radial_tree(
        self,
        node: _TreeNode,
        segments: list[str],
        labels: list[str],
        *,
        metadata: dict[str, dict[str, Any]],
        center_x: float,
        center_y: float,
    ) -> None:
        if node.left is not None and node.right is not None:
            segments.append(
                f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="2.8" fill="#f3efe6" stroke="#5d7268" stroke-width="1"></circle>'
            )
            for child in (node.left, node.right):
                segments.append(
                    f'<line x1="{node.x:.1f}" y1="{node.y:.1f}" x2="{child.x:.1f}" y2="{child.y:.1f}" stroke="#5d7268" stroke-width="1.9"></line>'
                )
                self._render_radial_tree(child, segments, labels, metadata=metadata, center_x=center_x, center_y=center_y)
            return

        species_label = node.name or ""
        meta = metadata.get(species_label, {})
        residue_fill, residue_text = self._residue_badge_colors(meta.get("match_class", "unknown"), species_label == "Human")
        dx = node.x - center_x
        dy = node.y - center_y
        mag = max((dx * dx + dy * dy) ** 0.5, 1.0)
        label_x = node.x + ((dx / mag) * 30.0)
        label_y = node.y + ((dy / mag) * 30.0)
        anchor = "start" if dx >= 0 else "end"
        badge_x = label_x + (8 if dx >= 0 else -52)
        identity_x = badge_x + (52 if dx >= 0 else -64)
        label_weight = "700" if species_label == "Human" else "500"
        label_fill = "#114f43" if species_label == "Human" else "#1e2b26"
        labels.append(
            f'<circle cx="{node.x:.1f}" cy="{node.y:.1f}" r="5.0" fill="{residue_fill}" stroke="#ffffff" stroke-width="1.4"></circle>'
            f'<text x="{label_x:.1f}" y="{label_y:.1f}" text-anchor="{anchor}" font-size="12" font-weight="{label_weight}" fill="{label_fill}">{escape(species_label)}</text>'
            f'<rect x="{badge_x:.1f}" y="{label_y - 12:.1f}" width="42" height="22" rx="11" fill="{residue_fill}"></rect>'
            f'<text x="{badge_x + 21:.1f}" y="{label_y + 3:.1f}" text-anchor="middle" font-size="11" font-weight="700" fill="{residue_text}">{escape(str(meta.get("site_residue", "?")))}</text>'
            f'<text x="{identity_x:.1f}" y="{label_y + 3:.1f}" text-anchor="{anchor}" font-size="10.5" fill="#5d7268">{float(meta.get("percent_identity", 0.0)):.1f}%</text>'
        )

    def _leaf_positions(self, root: _TreeNode) -> dict[str, float]:
        positions: dict[str, float] = {}
        self._collect_leaf_positions(root, positions)
        return positions

    def _collect_leaf_positions(self, node: _TreeNode, positions: dict[str, float]) -> None:
        if node.name is not None:
            positions[node.name] = node.y
            return
        if node.left is not None:
            self._collect_leaf_positions(node.left, positions)
        if node.right is not None:
            self._collect_leaf_positions(node.right, positions)

    @staticmethod
    def _residue_badge_colors(match_class: str, is_human: bool) -> tuple[str, str]:
        if is_human:
            return "#114f43", "#f7fbf8"
        if match_class == "match":
            return "#d6efe3", "#114f43"
        if match_class == "mismatch":
            return "#f5d8d0", "#8f3e28"
        return "#e9ece8", "#4f615a"

    @staticmethod
    def _identity_badge_colors(identity_pct: float, is_human: bool) -> tuple[str, str]:
        if is_human:
            return "#114f43", "#f7fbf8"
        if identity_pct >= 97.0:
            return "#d6efe3", "#114f43"
        if identity_pct >= 92.0:
            return "#dbe7f5", "#214a69"
        if identity_pct >= 85.0:
            return "#f4e7cf", "#8a5d16"
        return "#f2dce6", "#8f3557"

    @staticmethod
    def _tree_scale_bar(*, max_height: float, x: float, y: float, width: float) -> str:
        scale_value = max(round(max_height / 3.0, 2), 0.01)
        scale_width = min(max((scale_value / max_height) * width, 36.0), width)
        return (
            f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x + scale_width:.1f}" y2="{y:.1f}" stroke="#5d7268" stroke-width="2.4"></line>'
            f'<line x1="{x:.1f}" y1="{y - 5:.1f}" x2="{x:.1f}" y2="{y + 5:.1f}" stroke="#5d7268" stroke-width="2"></line>'
            f'<line x1="{x + scale_width:.1f}" y1="{y - 5:.1f}" x2="{x + scale_width:.1f}" y2="{y + 5:.1f}" stroke="#5d7268" stroke-width="2"></line>'
            f'<text x="{x + scale_width + 10:.1f}" y="{y + 4:.1f}" font-size="11" fill="#5d7268">branch length {scale_value:.2f}</text>'
        )

    def _newick_from_tree(self, node: _TreeNode) -> str:
        if node.name is not None:
            return self._sanitize_newick_label(node.name)
        left_length = max(node.height - (node.left.height if node.left else 0.0), 0.0)
        right_length = max(node.height - (node.right.height if node.right else 0.0), 0.0)
        left_text = self._newick_from_tree(node.left) if node.left else ""
        right_text = self._newick_from_tree(node.right) if node.right else ""
        return f"({left_text}:{left_length:.4f},{right_text}:{right_length:.4f})"

    @staticmethod
    def _sanitize_newick_label(label: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
