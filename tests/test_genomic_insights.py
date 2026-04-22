from __future__ import annotations

import requests
import unittest
from unittest.mock import patch

from app.services.genomic_insights import GenomicInsightService
from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import parse_missense_variant


class GenomicInsightTests(unittest.TestCase):
    def test_collect_builds_conservation_and_chromosome_panels(self) -> None:
        protein = ProteinRecord(
            accession="P00001",
            entry_name="TEST_HUMAN",
            gene_symbol="TEST",
            protein_name="Test protein",
            organism_name="Homo sapiens",
            sequence="MARRK",
            features=[],
            refseq_proteins=[],
            synonyms=[],
        )
        variant = parse_missense_variant("R3W")

        def fake_request(path: str, params=None):
            if path.startswith("/lookup/symbol/homo_sapiens/TEST"):
                return {
                    "id": "ENSG000001",
                    "display_name": "TEST",
                    "seq_region_name": "1",
                    "start": 100,
                    "end": 900,
                    "strand": 1,
                    "Transcript": [
                        {
                            "id": "ENST000001",
                            "start": 100,
                            "end": 900,
                            "strand": 1,
                            "is_canonical": 1,
                            "Translation": {"id": "ENSP000001", "length": 5},
                            "Exon": [
                                {"start": 100, "end": 180},
                                {"start": 400, "end": 500},
                                {"start": 800, "end": 900},
                            ],
                        }
                    ],
                }
            if path.startswith("/info/assembly/homo_sapiens/1"):
                return {"length": 1000}
            if path.startswith("/map/translation/ENSP000001/3..3"):
                return {"mappings": [{"mapped": {"seq_region_name": "1", "start": 410, "end": 412}}]}
            if path.startswith("/homology/symbol/human/TEST"):
                return {
                    "data": [
                        {
                            "id": "ENSG000001",
                            "homologies": [
                                {
                                    "source": {"align_seq": "MARRK"},
                                    "target": {
                                        "species": "pan_troglodytes",
                                        "align_seq": "MARRK",
                                        "perc_id": 99.1,
                                    },
                                },
                                {
                                    "source": {"align_seq": "MARRK"},
                                    "target": {
                                        "species": "macaca_mulatta",
                                        "align_seq": "MAQRK",
                                        "perc_id": 92.4,
                                    },
                                },
                            ],
                        }
                    ]
                }
            raise AssertionError(f"Unexpected path: {path}")

        service = GenomicInsightService(
            ensembl_base_url="https://rest.ensembl.org",
            ncbi_datasets_base_url="https://api.ncbi.nlm.nih.gov/datasets/v2",
            timeout_seconds=5.0,
            user_agent="test-agent",
        )

        with patch.object(service, "_request_json", side_effect=fake_request):
            bundle = service.collect(protein, variant)

        self.assertEqual(bundle["chromosome_location_panel"]["chromosome"], "1")
        self.assertEqual(bundle["chromosome_location_panel"]["exon_focus"]["label"], "Exon 2")
        self.assertEqual(bundle["conservation_panel"]["conserved_count"], 2)
        self.assertEqual(bundle["conservation_panel"]["dominant_residue"], "R")
        self.assertEqual(bundle["conservation_panel"]["mutant_match_count"], 0)
        self.assertEqual(bundle["phylogenetic_panel"]["closest_species"]["species_label"], "Chimpanzee")
        self.assertIn("<svg", bundle["phylogenetic_panel"]["tree_svg"])
        self.assertTrue(bundle["phylogenetic_panel"]["newick"].endswith(";"))
        self.assertEqual(bundle["phylogenetic_panel"]["matrix_rows"][0]["cells"][0]["class_name"], "matrix-self")

    def test_collect_uses_ncbi_fallback_for_chromosome_panel(self) -> None:
        protein = ProteinRecord(
            accession="P04637",
            entry_name="P53_HUMAN",
            gene_symbol="TP53",
            protein_name="Cellular tumor antigen p53",
            organism_name="Homo sapiens",
            sequence="MEEPQSDPSV",
            features=[],
            refseq_proteins=["NP_000537"],
            synonyms=[],
        )
        variant = parse_missense_variant("Q5Y")

        ensembl_service = GenomicInsightService(
            ensembl_base_url="https://rest.ensembl.org",
            ncbi_datasets_base_url="https://api.ncbi.nlm.nih.gov/datasets/v2",
            timeout_seconds=5.0,
            user_agent="test-agent",
        )

        ncbi_payload = {
            "reports": [
                {
                    "geneId": "7157",
                    "symbol": "TP53",
                    "transcripts": [
                        {
                            "accessionVersion": "NM_000546.6",
                            "ensemblTranscript": "ENST00000269305.9",
                            "type": "PROTEIN_CODING",
                            "cds": {"range": [{"begin": "1", "end": "30", "orientation": "plus"}]},
                            "protein": {
                                "accessionVersion": "NP_000537.3",
                                "length": 10,
                                "ensemblProtein": "ENSP00000269305.4",
                            },
                            "genomicLocations": [
                                {
                                    "genomicAccessionVersion": "NC_000017.11",
                                    "sequenceName": "Chromosome 17 Reference GRCh38.p14 Primary Assembly",
                                    "genomicRange": {"begin": "7668402", "end": "7687550", "orientation": "minus"},
                                    "exons": [
                                        {"begin": "7687376", "end": "7687550", "orientation": "minus", "order": 1},
                                        {"begin": "7676521", "end": "7676594", "orientation": "minus", "order": 2},
                                        {"begin": "7668402", "end": "7669690", "orientation": "minus", "order": 3},
                                    ],
                                }
                            ],
                        }
                    ],
                }
            ]
        }

        with patch.object(
            ensembl_service,
            "_request_json",
            side_effect=requests.RequestException("ensembl unavailable"),
        ), patch.object(
            ensembl_service,
            "_request_ncbi_json",
            return_value=ncbi_payload,
        ), patch.object(
            ensembl_service,
            "_build_conservation_panel",
            return_value=None,
        ):
            bundle = ensembl_service.collect(protein, variant)

        self.assertEqual(bundle["chromosome_location_panel"]["source"], "NCBI Datasets")
        self.assertEqual(bundle["chromosome_location_panel"]["chromosome"], "17")
        self.assertEqual(bundle["chromosome_location_panel"]["strand_label"], "reverse")
        self.assertTrue(any("NCBI Datasets fallback" in notice["body"] for notice in bundle["insight_notices"]))


if __name__ == "__main__":
    unittest.main()
