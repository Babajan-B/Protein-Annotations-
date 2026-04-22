from __future__ import annotations

import unittest

from app.services.variant_parser import VariantParseError, parse_missense_variant


class VariantParserTests(unittest.TestCase):
    def test_parses_short_notation(self) -> None:
        variant = parse_missense_variant("R248W")
        self.assertEqual(variant.wild_type, "R")
        self.assertEqual(variant.position, 248)
        self.assertEqual(variant.mutant, "W")
        self.assertEqual(variant.hgvs_protein, "p.Arg248Trp")

    def test_parses_hgvs_notation(self) -> None:
        variant = parse_missense_variant("p.Arg248Trp")
        self.assertEqual(variant.short_notation, "R248W")

    def test_rejects_identity_substitution(self) -> None:
        with self.assertRaises(VariantParseError):
            parse_missense_variant("R248R")

    def test_rejects_unsupported_format(self) -> None:
        with self.assertRaises(VariantParseError):
            parse_missense_variant("c.743G>A")


if __name__ == "__main__":
    unittest.main()
