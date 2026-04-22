from __future__ import annotations

import unittest
from unittest.mock import patch

from app.services.annotation_sources import AnnotationAggregator
from app.services.protein_lookup import ProteinRecord
from app.services.variant_parser import parse_missense_variant


class AnnotationAggregatorTests(unittest.TestCase):
    def test_collect_uses_uniprot_features_without_network(self) -> None:
        protein = ProteinRecord(
            accession="P04637",
            entry_name="P53_HUMAN",
            gene_symbol="TP53",
            protein_name="Cellular tumor antigen p53",
            organism_name="Homo sapiens",
            sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
            features=[
                {
                    "type": "DOMAIN",
                    "description": "Transactivation domain",
                    "location": {"start": {"value": 1}, "end": {"value": 20}},
                },
                {
                    "type": "MOD_RES",
                    "description": "Phosphoserine",
                    "location": {"start": {"value": 15}, "end": {"value": 15}},
                    "evidences": [{"evidenceCode": "ECO:0000269"}],
                },
            ],
            refseq_proteins=["NP_000537"],
            synonyms=[],
        )
        variant = parse_missense_variant("S15W")
        aggregator = AnnotationAggregator(
            interpro_base_url="https://example.org/interpro",
            ebi_proteins_base_url="https://example.org/proteins",
            ensembl_base_url="https://example.org/ensembl",
            ncbi_base_url="https://example.org/ncbi",
            mavedb_base_url="https://example.org/mavedb",
            timeout_seconds=0.1,
            user_agent="test-agent",
        )

        with patch.object(aggregator, "_request_json", return_value=None):
            annotations = aggregator.collect(protein, variant)

        self.assertEqual(len(annotations["interval_features"]), 1)
        self.assertEqual(len(annotations["site_features"]), 1)
        self.assertEqual(annotations["site_features"][0]["position"], 15)
        self.assertTrue(any(row["label"] == "Domains and motifs" for row in annotations["feature_counts"]))
        self.assertTrue(any(status["source"] == "InterPro API" for status in annotations["source_status"]))


if __name__ == "__main__":
    unittest.main()
