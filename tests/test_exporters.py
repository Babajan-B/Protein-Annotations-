from __future__ import annotations

import unittest

from app.services.analysis import MetricDelta
from app.services.exporters import report_to_csv, serialize_report, structure_annotations_to_csv
from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import parse_missense_variant


class ExporterTests(unittest.TestCase):
    def test_json_serialization_and_csv_export(self) -> None:
        report = {
            "protein": ProteinRecord(
                accession="P00001",
                entry_name="TEST_HUMAN",
                gene_symbol="TEST",
                protein_name="Test protein",
                organism_name="Homo sapiens",
                sequence="MARRK",
                features=[],
                refseq_proteins=[],
                synonyms=[],
            ),
            "variant": parse_missense_variant("R3W"),
            "metric_deltas": [MetricDelta(label="GRAVY", wild_type=1.0, mutant=2.0, delta=1.0)],
            "interval_features": [],
            "site_features": [],
            "evidence_rows": [],
            "predictor_tables": [],
            "residue_tracks": [],
            "resource_cards": [],
            "source_status": [],
            "accessibility_panel": None,
        }

        serialized = serialize_report(report)
        csv_output = report_to_csv(report)

        self.assertEqual(serialized["protein"]["accession"], "P00001")
        self.assertNotIn("predictor_tables", serialized)
        self.assertNotIn("resource_cards", serialized)
        self.assertNotIn("source_status", serialized)
        self.assertNotIn("accessibility_panel", serialized)
        self.assertIn("GRAVY", csv_output)
        self.assertIn("TEST_HUMAN", csv_output)

    def test_structure_annotations_csv_export(self) -> None:
        report = {
            "secondary_structure_panel": {
                "annotation_blocks": [
                    {
                        "start": 1,
                        "end": 2,
                        "residues": [
                            {
                                "position": 1,
                                "aa": "M",
                                "state": "Coil",
                                "state_code": "C",
                                "confidence": 0.41,
                                "coil_probability": 0.41,
                                "helix_probability": 0.20,
                                "strand_probability": 0.39,
                                "is_variant_site": False,
                            },
                            {
                                "position": 2,
                                "aa": "A",
                                "state": "Alpha helix",
                                "state_code": "H",
                                "confidence": 0.82,
                                "coil_probability": 0.06,
                                "helix_probability": 0.82,
                                "strand_probability": 0.12,
                                "is_variant_site": True,
                            },
                        ],
                    }
                ]
            }
        }

        csv_output = structure_annotations_to_csv(report)

        self.assertIn("position,amino_acid,state,state_code", csv_output)
        self.assertIn("2,A,Alpha helix,H,0.82,0.06,0.82,0.12,yes", csv_output)


if __name__ == "__main__":
    unittest.main()
