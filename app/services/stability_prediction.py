from __future__ import annotations

from typing import Any

# BLOSUM62 substitution matrix (selected rows for the 20 standard amino acids)
# Values from the standard BLOSUM62 matrix
BLOSUM62: dict[str, dict[str, int]] = {
    "A": {"A": 4,"R":-1,"N":-2,"D":-2,"C": 0,"Q":-1,"E":-1,"G": 0,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S": 1,"T": 0,"W":-3,"Y":-2,"V": 0},
    "R": {"A":-1,"R": 5,"N": 0,"D":-2,"C":-3,"Q": 1,"E": 0,"G":-2,"H": 0,"I":-3,"L":-2,"K": 2,"M":-1,"F":-3,"P":-2,"S":-1,"T":-1,"W":-3,"Y":-2,"V":-3},
    "N": {"A":-2,"R": 0,"N": 6,"D": 1,"C":-3,"Q": 0,"E": 0,"G": 0,"H": 1,"I":-3,"L":-3,"K": 0,"M":-2,"F":-3,"P":-2,"S": 1,"T": 0,"W":-4,"Y":-2,"V":-3},
    "D": {"A":-2,"R":-2,"N": 1,"D": 6,"C":-3,"Q": 0,"E": 2,"G":-1,"H":-1,"I":-3,"L":-4,"K":-1,"M":-3,"F":-3,"P":-1,"S": 0,"T":-1,"W":-4,"Y":-3,"V":-3},
    "C": {"A": 0,"R":-3,"N":-3,"D":-3,"C": 9,"Q":-3,"E":-4,"G":-3,"H":-3,"I":-1,"L":-1,"K":-3,"M":-1,"F":-2,"P":-3,"S":-1,"T":-1,"W":-2,"Y":-2,"V":-1},
    "Q": {"A":-1,"R": 1,"N": 0,"D": 0,"C":-3,"Q": 5,"E": 2,"G":-2,"H": 0,"I":-3,"L":-2,"K": 1,"M": 0,"F":-3,"P":-1,"S": 0,"T":-1,"W":-2,"Y":-1,"V":-2},
    "E": {"A":-1,"R": 0,"N": 0,"D": 2,"C":-4,"Q": 2,"E": 5,"G":-2,"H": 0,"I":-3,"L":-3,"K": 1,"M":-2,"F":-3,"P":-1,"S": 0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "G": {"A": 0,"R":-2,"N": 0,"D":-1,"C":-3,"Q":-2,"E":-2,"G": 6,"H":-2,"I":-4,"L":-4,"K":-2,"M":-3,"F":-3,"P":-2,"S": 0,"T":-2,"W":-2,"Y":-3,"V":-3},
    "H": {"A":-2,"R": 0,"N": 1,"D":-1,"C":-3,"Q": 0,"E": 0,"G":-2,"H": 8,"I":-3,"L":-3,"K":-1,"M":-2,"F":-1,"P":-2,"S":-1,"T":-2,"W":-2,"Y": 2,"V":-3},
    "I": {"A":-1,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-3,"E":-3,"G":-4,"H":-3,"I": 4,"L": 2,"K":-3,"M": 1,"F": 0,"P":-3,"S":-2,"T":-1,"W":-3,"Y":-1,"V": 3},
    "L": {"A":-1,"R":-2,"N":-3,"D":-4,"C":-1,"Q":-2,"E":-3,"G":-4,"H":-3,"I": 2,"L": 4,"K":-2,"M": 2,"F": 0,"P":-3,"S":-2,"T":-1,"W":-2,"Y":-1,"V": 1},
    "K": {"A":-1,"R": 2,"N": 0,"D":-1,"C":-3,"Q": 1,"E": 1,"G":-2,"H":-1,"I":-3,"L":-2,"K": 5,"M":-1,"F":-3,"P":-1,"S": 0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "M": {"A":-1,"R":-1,"N":-2,"D":-3,"C":-1,"Q": 0,"E":-2,"G":-3,"H":-2,"I": 1,"L": 2,"K":-1,"M": 5,"F": 0,"P":-2,"S":-1,"T":-1,"W":-1,"Y":-1,"V": 1},
    "F": {"A":-2,"R":-3,"N":-3,"D":-3,"C":-2,"Q":-3,"E":-3,"G":-3,"H":-1,"I": 0,"L": 0,"K":-3,"M": 0,"F": 6,"P":-4,"S":-2,"T":-2,"W": 1,"Y": 3,"V":-1},
    "P": {"A":-1,"R":-2,"N":-2,"D":-1,"C":-3,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-3,"L":-3,"K":-1,"M":-2,"F":-4,"P": 7,"S":-1,"T":-1,"W":-4,"Y":-3,"V":-2},
    "S": {"A": 1,"R":-1,"N": 1,"D": 0,"C":-1,"Q": 0,"E": 0,"G": 0,"H":-1,"I":-2,"L":-2,"K": 0,"M":-1,"F":-2,"P":-1,"S": 4,"T": 1,"W":-3,"Y":-2,"V":-2},
    "T": {"A": 0,"R":-1,"N": 0,"D":-1,"C":-1,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S": 1,"T": 5,"W":-2,"Y":-2,"V": 0},
    "W": {"A":-3,"R":-3,"N":-4,"D":-4,"C":-2,"Q":-2,"E":-3,"G":-2,"H":-2,"I":-3,"L":-2,"K":-3,"M":-1,"F": 1,"P":-4,"S":-3,"T":-2,"W":11,"Y": 2,"V":-3},
    "Y": {"A":-2,"R":-2,"N":-2,"D":-3,"C":-2,"Q":-1,"E":-2,"G":-3,"H": 2,"I":-1,"L":-1,"K":-2,"M":-1,"F": 3,"P":-3,"S":-2,"T":-2,"W": 2,"Y": 7,"V":-1},
    "V": {"A": 0,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-2,"E":-2,"G":-3,"H":-3,"I": 3,"L": 1,"K":-2,"M": 1,"F":-1,"P":-2,"S":-2,"T": 0,"W":-3,"Y":-1,"V": 4},
}

# Hydrophobicity scale (Kyte-Doolittle, same as analysis.py)
_HYDROPATHY: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Residue volume (Å³, approximate)
_VOLUME: dict[str, float] = {
    "G": 60, "A": 89, "S": 96, "C": 108, "T": 116, "P": 112, "V": 140,
    "D": 111, "N": 114, "I": 167, "L": 167, "E": 138, "Q": 144, "M": 163,
    "H": 153, "F": 190, "R": 173, "Y": 193, "W": 228, "K": 168,
}

_CHARGE: dict[str, str] = {
    "R": "positive", "K": "positive", "H": "positive",
    "D": "negative", "E": "negative",
}

_POLARITY: dict[str, str] = {
    "R": "polar", "N": "polar", "D": "polar", "C": "polar", "Q": "polar",
    "E": "polar", "H": "polar", "K": "polar", "S": "polar", "T": "polar", "Y": "polar",
}

_AROMATIC = frozenset("FWYH")


class StabilityPredictor:
    """Sequence-based mutation stability estimator.

    Uses BLOSUM62 evolutionary compatibility score and physicochemical
    property differences to estimate mutation impact on protein stability.
    ΔΔG is an empirical approximation, not a thermodynamic calculation.
    """

    def collect(
        self,
        *,
        wild_type: str,
        mutant: str,
        position: int,
    ) -> dict[str, Any]:
        wt = wild_type.upper()
        mt = mutant.upper()

        if wt not in BLOSUM62 or mt not in BLOSUM62:
            return {"stability_panel": None}

        blosum_score = BLOSUM62[wt].get(mt, 0)
        blosum_self = BLOSUM62[wt].get(wt, 1)

        # Normalise BLOSUM62 score against the self-substitution score
        blosum_norm = blosum_score / max(abs(blosum_self), 1)

        # Hydrophobicity delta
        hydro_delta = _HYDROPATHY.get(mt, 0.0) - _HYDROPATHY.get(wt, 0.0)

        # Volume delta (normalised to ±1 range via max possible ~168 Å³)
        vol_wt = _VOLUME.get(wt, 120.0)
        vol_mt = _VOLUME.get(mt, 120.0)
        vol_delta_norm = (vol_mt - vol_wt) / 168.0

        # Empirical ΔΔG estimate (kcal/mol)
        # Formula: penalise BLOSUM mismatches + hydrophobicity disruption + volume change
        ddg = round(
            -0.8 * blosum_norm       # evolutionary incompatibility → destabilising
            + 0.3 * abs(hydro_delta)  # large hydrophobicity shift → destabilising
            + 0.2 * abs(vol_delta_norm),  # large size change → destabilising
            2,
        )

        if ddg < -0.5:
            ddg_class = "stabilizing"
            ddg_label = "Stabilizing"
            ddg_class_name = "ddg-stabilizing"
        elif ddg > 0.5:
            ddg_class = "destabilizing"
            ddg_label = "Destabilizing"
            ddg_class_name = "ddg-destabilizing"
        else:
            ddg_class = "neutral"
            ddg_label = "Neutral"
            ddg_class_name = "ddg-neutral"

        # BLOSUM62 interpretation
        if blosum_score >= 2:
            blosum_interpretation = "Conservative substitution (evolutionary compatible)"
        elif blosum_score >= 0:
            blosum_interpretation = "Moderately conservative substitution"
        elif blosum_score >= -2:
            blosum_interpretation = "Semi-conservative substitution (mild penalty)"
        else:
            blosum_interpretation = "Non-conservative substitution (strong evolutionary penalty)"

        # Physicochemical property changes
        property_changes = _property_changes(wt, mt)

        # SIFT-proxy: fraction of BLOSUM62 row entries that are ≥ blosum_score
        # (higher = more positions tolerate this kind of change)
        row_values = list(BLOSUM62[wt].values())
        sift_proxy = round(sum(1 for v in row_values if v >= blosum_score) / len(row_values), 3)
        sift_tolerance = "Tolerated" if sift_proxy >= 0.5 else "Intolerant"

        summary = (
            f"The {wt}{position}{mt} substitution has a BLOSUM62 score of {blosum_score} "
            f"({blosum_interpretation.lower()}) and an estimated ΔΔG of {ddg:+.2f} kcal/mol "
            f"({ddg_label.lower()}). "
            f"SIFT-proxy tolerance score: {sift_proxy} ({sift_tolerance})."
        )

        return {
            "stability_panel": {
                "ddg_estimate": ddg,
                "ddg_class": ddg_class,
                "ddg_label": ddg_label,
                "ddg_class_name": ddg_class_name,
                "blosum62_score": blosum_score,
                "blosum62_interpretation": blosum_interpretation,
                "hydro_delta": round(hydro_delta, 2),
                "vol_delta": round(vol_mt - vol_wt, 1),
                "sift_proxy": sift_proxy,
                "sift_tolerance": sift_tolerance,
                "property_changes": property_changes,
                "summary": summary,
                "disclaimer": (
                    "ΔΔG is an empirical approximation from sequence-based features "
                    "(BLOSUM62 + hydrophobicity + volume), not a physics-based calculation."
                ),
            }
        }


def _property_changes(wt: str, mt: str) -> list[dict[str, str]]:
    changes = []

    wt_charge = _CHARGE.get(wt, "neutral")
    mt_charge = _CHARGE.get(mt, "neutral")
    if wt_charge != mt_charge:
        changes.append({
            "property": "Charge",
            "from": wt_charge,
            "to": mt_charge,
            "impact": "significant",
        })

    wt_polar = _POLARITY.get(wt, "non-polar")
    mt_polar = _POLARITY.get(mt, "non-polar")
    if wt_polar != mt_polar:
        changes.append({
            "property": "Polarity",
            "from": wt_polar,
            "to": mt_polar,
            "impact": "moderate",
        })

    wt_arom = wt in _AROMATIC
    mt_arom = mt in _AROMATIC
    if wt_arom != mt_arom:
        changes.append({
            "property": "Aromaticity",
            "from": "aromatic" if wt_arom else "non-aromatic",
            "to": "aromatic" if mt_arom else "non-aromatic",
            "impact": "moderate",
        })

    vol_diff = abs(_VOLUME.get(mt, 120.0) - _VOLUME.get(wt, 120.0))
    if vol_diff > 40:
        changes.append({
            "property": "Volume",
            "from": f"{_VOLUME.get(wt, 120.0):.0f} Å³",
            "to": f"{_VOLUME.get(mt, 120.0):.0f} Å³",
            "impact": "significant" if vol_diff > 80 else "moderate",
        })

    if not changes:
        changes.append({
            "property": "All compared properties",
            "from": "unchanged",
            "to": "unchanged",
            "impact": "minimal",
        })

    return changes
