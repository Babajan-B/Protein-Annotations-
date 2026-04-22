from __future__ import annotations

import unittest

from app.services.analysis import SequenceAnalysisService
from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import parse_missense_variant


class AnalysisServiceTests(unittest.TestCase):
    def test_build_report_with_fallback_metrics(self) -> None:
        protein = ProteinRecord(
            accession="P00001",
            entry_name="TEST_HUMAN",
            gene_symbol="TEST",
            protein_name="Test protein",
            organism_name="Homo sapiens",
            sequence="MARRKSTVVVVVVVVV",
            features=[],
            refseq_proteins=["NP_000001"],
            synonyms=["Test protein alt"],
        )
        variant = parse_missense_variant("R3W")
        annotations = {
            "interval_features": [
                {
                    "source": "UniProt",
                    "type": "DOMAIN",
                    "category": "domains",
                    "label": "DNA-binding",
                    "description": "DNA-binding",
                    "start": 2,
                    "end": 6,
                    "evidence": "",
                    "variant_hit": True,
                }
            ],
            "site_features": [
                {
                    "source": "UniProt",
                    "type": "MOD_RES",
                    "category": "ptm",
                    "label": "Phosphoserine",
                    "description": "Phosphoserine",
                    "position": 3,
                    "start": 3,
                    "end": 3,
                    "evidence": "ECO:0000269",
                    "variant_hit": True,
                }
            ],
            "evidence_rows": [
                {
                    "source": "ClinVar",
                    "kind": "clinvar",
                    "label": "Example record",
                    "position": 3,
                    "summary": "Pathogenic",
                    "link_hint": "RCV000000",
                    "score": "",
                }
            ],
            "architecture_rows": [],
            "feature_counts": [{"label": "Domains and motifs", "count": 1}],
            "evidence_plot": [{"position": 3, "left_pct": 18.75, "hit_variant_site": True}],
            "source_status": [
                {"source": "ClinVar", "status": "ok", "summary": "Loaded example"},
                {"source": "DMVFL-RSA", "status": "error", "summary": "Timed out"},
            ],
            "resource_cards": [
                {
                    "key": "clinvar",
                    "label": "ClinVar",
                    "category": "evidence",
                    "summary": "Clinical evidence",
                    "command_env": "NA",
                    "configured": "configured",
                    "docs_url": "https://example.test/clinvar",
                    "resource_hint": "Loaded example",
                },
                {
                    "key": "dmvfl_rsa",
                    "label": "DMVFL-RSA",
                    "category": "structure",
                    "summary": "Accessibility",
                    "command_env": "NA",
                    "configured": "configured",
                    "docs_url": "https://example.test/dmvfl",
                    "resource_hint": "Timed out",
                },
            ],
            "predictor_tables": [
                {
                    "title": "PSIPRED residue assignments",
                    "columns": [
                        "Position",
                        "AA",
                        "State",
                        "State code",
                        "Confidence",
                        "Coil probability",
                        "Helix probability",
                        "Strand probability",
                    ],
                    "rows": [
                        {
                            "Position": 1,
                            "AA": "M",
                            "State": "Coil",
                            "State code": "C",
                            "Confidence": 0.66,
                            "Coil probability": 0.66,
                            "Helix probability": 0.12,
                            "Strand probability": 0.22,
                        },
                        {
                            "Position": 2,
                            "AA": "A",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.88,
                            "Coil probability": 0.02,
                            "Helix probability": 0.88,
                            "Strand probability": 0.10,
                        },
                        {
                            "Position": 3,
                            "AA": "R",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.91,
                            "Coil probability": 0.03,
                            "Helix probability": 0.91,
                            "Strand probability": 0.06,
                        },
                    ],
                },
                {
                    "title": "DMVFL-RSA residue accessibility",
                    "columns": ["Position", "AA", "RSA", "ASA", "Exposure"],
                    "rows": [
                        {"Position": 1, "AA": "M", "RSA": 0.11, "ASA": 20.0, "Exposure": "Intermediate"},
                        {"Position": 2, "AA": "A", "RSA": 0.34, "ASA": 39.1, "Exposure": "Exposed"},
                        {"Position": 3, "AA": "R", "RSA": 0.08, "ASA": 18.0, "Exposure": "Buried"},
                    ],
                }
            ],
            "notices": [],
        }

        service = SequenceAnalysisService(
            [
                {"key": "resolver", "label": "Resolver", "status": "implemented", "summary": "ok"},
                {"key": "evidence", "label": "Evidence", "status": "partial", "summary": "ok"},
            ]
        )
        report = service.build_report(protein=protein, variant=variant, annotations=annotations)

        self.assertEqual(report["variant"].short_notation, "R3W")
        self.assertEqual(report["context"]["wild_window"][2], "R")
        self.assertEqual(report["interval_features"][0]["label"], "DNA-binding")
        self.assertEqual(report["summary_cards"][1]["label"], "Feature overlap")
        self.assertTrue(any(item["label"] == "GRAVY" for item in report["metric_plot"]))
        self.assertEqual(report["accessibility_panel"]["variant_row"]["Exposure"], "Buried")
        self.assertEqual(report["secondary_structure_panel"]["variant_row"]["state"], "Alpha helix")
        self.assertEqual(report["secondary_structure_panel"]["neighborhood_rows"][2]["position"], 3)
        self.assertEqual(report["source_status"], [{"source": "ClinVar", "status": "ok", "summary": "Loaded example"}])
        self.assertEqual(report["resource_cards"][0]["label"], "ClinVar")
        self.assertEqual(report["supplementary_predictor_tables"], [])


if __name__ == "__main__":
    unittest.main()
