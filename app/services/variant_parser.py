from __future__ import annotations

from dataclasses import dataclass
import re


ONE_TO_THREE = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}

THREE_TO_ONE = {value.upper(): key for key, value in ONE_TO_THREE.items()}

SHORT_VARIANT_RE = re.compile(
    r"^\s*(?P<wt>[ACDEFGHIKLMNPQRSTVWY])\s*(?P<pos>\d+)\s*(?P<mut>[ACDEFGHIKLMNPQRSTVWY])\s*$",
    re.IGNORECASE,
)
HGVS_VARIANT_RE = re.compile(
    r"^\s*p\.\s*(?P<wt>[A-Za-z]{3})\s*(?P<pos>\d+)\s*(?P<mut>[A-Za-z]{3})\s*$",
    re.IGNORECASE,
)


class VariantParseError(ValueError):
    """Raised when the user submits an unsupported protein variant."""


@dataclass(frozen=True)
class MissenseVariant:
    original_text: str
    wild_type: str
    position: int
    mutant: str

    @property
    def short_notation(self) -> str:
        return f"{self.wild_type}{self.position}{self.mutant}"

    @property
    def hgvs_protein(self) -> str:
        return f"p.{ONE_TO_THREE[self.wild_type]}{self.position}{ONE_TO_THREE[self.mutant]}"


def parse_missense_variant(text: str) -> MissenseVariant:
    candidate = (text or "").strip()
    if not candidate:
        raise VariantParseError("Variant is required. Use a missense substitution such as R248W.")

    short_match = SHORT_VARIANT_RE.match(candidate)
    if short_match:
        wild_type = short_match.group("wt").upper()
        mutant = short_match.group("mut").upper()
        position = int(short_match.group("pos"))
        return _build_variant(candidate, wild_type, position, mutant)

    hgvs_match = HGVS_VARIANT_RE.match(candidate)
    if hgvs_match:
        wild_type = _three_letter_to_one(hgvs_match.group("wt"))
        mutant = _three_letter_to_one(hgvs_match.group("mut"))
        position = int(hgvs_match.group("pos"))
        return _build_variant(candidate, wild_type, position, mutant)

    raise VariantParseError(
        "Unsupported variant format. Use either short notation like R248W or HGVS protein notation like p.Arg248Trp."
    )


def _build_variant(original_text: str, wild_type: str, position: int, mutant: str) -> MissenseVariant:
    if position <= 0:
        raise VariantParseError("Residue positions must be positive integers.")
    if wild_type == mutant:
        raise VariantParseError("The mutant residue must differ from the wild-type residue for missense analysis.")
    return MissenseVariant(
        original_text=original_text,
        wild_type=wild_type,
        position=position,
        mutant=mutant,
    )


def _three_letter_to_one(token: str) -> str:
    key = token.strip().upper()
    if key not in THREE_TO_ONE:
        raise VariantParseError(f"Unsupported amino-acid code '{token}'.")
    return THREE_TO_ONE[key]
