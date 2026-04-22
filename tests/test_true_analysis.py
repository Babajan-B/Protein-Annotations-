from __future__ import annotations

import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from app.services.protein_lookup import ProteinRecord
from app.services.tool_registry import build_resource_cards
from app.services.true_analysis import TrueAnalysisOrchestrator
from app.services.variant_parser import parse_missense_variant


class TrueAnalysisTests(unittest.TestCase):
    def test_resource_cards_reflect_env_configuration(self) -> None:
        with patch.dict(os.environ, {"IUPRED3_COMMAND_TEMPLATE": "python wrapper.py"}, clear=False):
            cards = build_resource_cards()
        iupred = next(card for card in cards if card["key"] == "iupred3")
        self.assertEqual(iupred["configured"], "configured")

    def test_orchestrator_runs_configured_wrapper(self) -> None:
        protein = ProteinRecord(
            accession="P00001",
            entry_name="TEST_HUMAN",
            gene_symbol="TEST",
            protein_name="Test protein",
            organism_name="Homo sapiens",
            sequence="MARRKSTVVVVVVVVV",
            features=[],
            refseq_proteins=[],
            synonyms=[],
        )
        variant = parse_missense_variant("R3W")
        script = Path(__file__).parent / "fixtures" / "mock_predictor.py"
        command = (
            f"python3 {script} --input-fasta {{input_fasta}} --output-json {{output_json}} "
            "--variant-position {variant_position}"
        )

        def fake_resolve_command_template(spec):
            if spec.command_env == "IUPRED3_COMMAND_TEMPLATE":
                return command
            return ""

        with tempfile.TemporaryDirectory() as tmpdir, patch(
            "app.services.true_analysis.resolve_command_template",
            side_effect=fake_resolve_command_template,
        ):
            orchestrator = TrueAnalysisOrchestrator(workdir=tmpdir, timeout_seconds=5.0)
            bundle = orchestrator.collect(protein, variant)

        self.assertEqual(bundle["interval_features"][0]["label"], "Mock domain")
        self.assertEqual(bundle["site_features"][0]["position"], 3)
        self.assertEqual(bundle["residue_tracks"][0]["label"], "Mock score")
        statuses = {row["source"]: row["status"] for row in bundle["source_status"]}
        self.assertEqual(statuses["IUPred3"], "ok")


if __name__ == "__main__":
    unittest.main()
