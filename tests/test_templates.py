from __future__ import annotations

import unittest

from flask import render_template

from app import create_app
from app.services.analysis import SequenceAnalysisService
from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import parse_missense_variant


class TemplateRenderTests(unittest.TestCase):
    def test_result_template_renders_residue_tracks(self) -> None:
        app = create_app()
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
            "interval_features": [],
            "site_features": [],
            "evidence_rows": [],
            "residue_tracks": [
                {
                    "key": "psipred_confidence",
                    "label": "PSIPRED confidence",
                    "source": "PSIPRED",
                    "values": [
                        {"position": 1, "value": 0.4, "height_pct": 40.0, "is_variant_site": False},
                        {"position": 2, "value": 0.8, "height_pct": 80.0, "is_variant_site": False},
                        {"position": 3, "value": 0.9, "height_pct": 90.0, "is_variant_site": True},
                    ],
                }
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
                            "Confidence": 0.41,
                            "Coil probability": 0.41,
                            "Helix probability": 0.20,
                            "Strand probability": 0.39,
                        },
                        {
                            "Position": 2,
                            "AA": "A",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.82,
                            "Coil probability": 0.06,
                            "Helix probability": 0.82,
                            "Strand probability": 0.12,
                        },
                        {
                            "Position": 3,
                            "AA": "R",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.91,
                            "Coil probability": 0.02,
                            "Helix probability": 0.91,
                            "Strand probability": 0.07,
                        },
                    ],
                }
            ],
            "resource_cards": [],
            "architecture_rows": [],
            "feature_counts": [],
            "evidence_plot": [],
            "source_status": [{"source": "PSIPRED", "status": "ok", "summary": "Loaded example"}],
            "notices": [],
        }
        service = SequenceAnalysisService(
            [
                {"key": "resolver", "label": "Resolver", "status": "implemented", "summary": "ok"},
                {"key": "structure", "label": "Structure", "status": "implemented", "summary": "ok"},
            ]
        )
        report = service.build_report(protein=protein, variant=variant, annotations=annotations)
        report["conservation_panel"] = {
            "conservation_pct": 87.5,
            "conserved_count": 7,
            "comparable_count": 8,
            "dominant_residue": "R",
            "human_residue": "R",
            "mutant_residue": "W",
            "mutant_match_count": 0,
            "variant_position": 3,
            "species_rows": [
                {"species_label": "Human", "site_residue": "R", "percent_identity": 100.0, "match_class": "match"}
            ],
        }
        report["phylogenetic_panel"] = {
            "method": "Whole-gene phylogeny approximated from full-length orthologue protein sequence identity using UPGMA clustering.",
            "species_count": 3,
            "closest_species": {"species_label": "Chimpanzee", "percent_identity": 99.1},
            "farthest_species": {"species_label": "Rhesus macaque", "percent_identity": 92.4},
            "matrix_headers": ["Human", "Chimpanzee", "Rhesus macaque"],
            "matrix_rows": [
                {
                    "species_label": "Human",
                    "cells": [
                        {"value": 100.0, "class_name": "matrix-self"},
                        {"value": 99.1, "class_name": "matrix-vhigh"},
                        {"value": 92.4, "class_name": "matrix-high"},
                    ],
                }
            ],
            "tree_svg": "<svg></svg>",
            "tree_views": [
                {"key": "publication", "label": "Publication", "summary": "Detailed rectangular tree.", "svg": "<svg></svg>"},
                {"key": "compact", "label": "Compact", "summary": "Dense rectangular tree.", "svg": "<svg></svg>"},
                {"key": "diagonal", "label": "Diagonal", "summary": "Slanted cladogram.", "svg": "<svg></svg>"},
                {"key": "radial", "label": "Radial", "summary": "Circular whole-gene tree.", "svg": "<svg></svg>"},
            ],
            "newick": "(Human:0.0000,(Chimpanzee:0.0100,Rhesus_macaque:0.0750):0.0200);",
        }
        report["chromosome_location_panel"] = {
            "chromosome": "1",
            "chromosome_length_label": "1,000 bp",
            "gene_symbol": "TEST",
            "gene_span_label": "100-900",
            "variant_span_label": "1:410-412",
            "strand_label": "forward",
            "transcript_id": "ENST000001",
            "gene_left_pct": 10.0,
            "gene_width_pct": 80.0,
            "variant_left_pct": 41.0,
            "gene_track_variant_pct": 38.0,
            "gene_track_variant_width_pct": 2.0,
            "exon_rows": [
                {"number": 1, "left_pct": 0.0, "width_pct": 24.0, "contains_variant": False},
                {"number": 2, "left_pct": 34.0, "width_pct": 28.0, "contains_variant": True},
            ],
            "exon_focus": {
                "label": "Exon 2",
                "start": 400,
                "end": 500,
                "length": 101,
                "visible_start": 400,
                "visible_end": 500,
                "left_truncated": False,
                "right_truncated": False,
                "variant_left_pct": 10.0,
                "variant_width_pct": 4.0,
            },
        }
        report["insight_notices"] = []

        with app.test_request_context("/report"):
            html = render_template("result.html", report=report, query_args={"query": "TEST", "variant": "R3W"})

        self.assertIn("Secondary structure annotation", html)
        self.assertIn("Structural annotations", html)
        self.assertIn("Variant neighborhood", html)
        self.assertIn("Primate residue conservation at the variant site", html)
        self.assertIn("Whole-gene phylogenetic analysis", html)
        self.assertIn("Phylogeny Newick", html)
        self.assertIn("Publication", html)
        self.assertIn("Compact", html)
        self.assertIn("Diagonal", html)
        self.assertIn("Radial", html)
        self.assertIn("Chromosome and transcript location", html)

    def test_structure_annotations_template_renders_annotation_blocks(self) -> None:
        app = create_app()
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
            "interval_features": [],
            "site_features": [],
            "evidence_rows": [],
            "residue_tracks": [],
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
                            "Confidence": 0.41,
                            "Coil probability": 0.41,
                            "Helix probability": 0.20,
                            "Strand probability": 0.39,
                        },
                        {
                            "Position": 2,
                            "AA": "A",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.82,
                            "Coil probability": 0.06,
                            "Helix probability": 0.82,
                            "Strand probability": 0.12,
                        },
                        {
                            "Position": 3,
                            "AA": "R",
                            "State": "Alpha helix",
                            "State code": "H",
                            "Confidence": 0.91,
                            "Coil probability": 0.02,
                            "Helix probability": 0.91,
                            "Strand probability": 0.07,
                        },
                    ],
                }
            ],
            "resource_cards": [],
            "architecture_rows": [],
            "feature_counts": [],
            "evidence_plot": [],
            "source_status": [{"source": "PSIPRED", "status": "ok", "summary": "Loaded example"}],
            "notices": [],
        }
        service = SequenceAnalysisService(
            [
                {"key": "resolver", "label": "Resolver", "status": "implemented", "summary": "ok"},
                {"key": "structure", "label": "Structure", "status": "implemented", "summary": "ok"},
            ]
        )
        report = service.build_report(protein=protein, variant=variant, annotations=annotations)

        with app.test_request_context("/structure-annotations"):
            html = render_template(
                "structure_annotations.html",
                report=report,
                query_args={"query": "TEST", "variant": "R3W"},
            )

        self.assertIn("Residue annotation blocks", html)
        self.assertIn("Alpha helix", html)
        self.assertIn("Position 3 carries amino acid", html)


if __name__ == "__main__":
    unittest.main()
