"""Variant enumeration and protein-change parsing."""

from __future__ import annotations

from dataclasses import dataclass
import re

from .constants import AMINO_ACIDS, MAPT_441_SEQUENCE

AA3_TO_AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}


@dataclass(frozen=True)
class MissenseVariant:
    position: int
    wt_aa: str
    mut_aa: str

    @property
    def variant_id(self) -> str:
        return f"{self.wt_aa}{self.position}{self.mut_aa}"

    @property
    def protein_change(self) -> str:
        return f"p.{self.wt_aa}{self.position}{self.mut_aa}"

    def as_row(self) -> dict[str, str | int]:
        return {
            "variant_id": self.variant_id,
            "protein_change": self.protein_change,
            "position": self.position,
            "wt_aa": self.wt_aa,
            "mut_aa": self.mut_aa,
        }


def validate_reference_sequence(sequence: str = MAPT_441_SEQUENCE) -> None:
    if len(sequence) != 441:
        raise ValueError(f"Expected 441-aa MAPT sequence, got {len(sequence)}")
    invalid = sorted(set(sequence) - set(AMINO_ACIDS))
    if invalid:
        raise ValueError(f"Invalid amino-acid symbols in reference sequence: {invalid}")


def generate_missense_variants(sequence: str = MAPT_441_SEQUENCE) -> list[MissenseVariant]:
    validate_reference_sequence(sequence)
    variants: list[MissenseVariant] = []
    for zero_based, wt_aa in enumerate(sequence):
        position = zero_based + 1
        for mut_aa in AMINO_ACIDS:
            if mut_aa != wt_aa:
                variants.append(MissenseVariant(position=position, wt_aa=wt_aa, mut_aa=mut_aa))
    return variants


def parse_protein_change(value: str) -> MissenseVariant | None:
    """Parse common protein change forms such as p.R406W or p.Arg406Trp."""
    if not value:
        return None

    text = value.strip()
    text = text.replace("(", "").replace(")", "")
    if ":" in text:
        text = text.split(":")[-1]
    text = text.removeprefix("NP_005901.2")
    text = text.removeprefix("p.")

    one_letter = re.search(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{1,4})([ACDEFGHIKLMNPQRSTVWY])\b", text)
    if one_letter:
        wt, pos, mut = one_letter.groups()
        return MissenseVariant(position=int(pos), wt_aa=wt, mut_aa=mut)

    three_letter = re.search(r"\b([A-Z][a-z]{2})(\d{1,4})([A-Z][a-z]{2})\b", text)
    if three_letter:
        wt3, pos, mut3 = three_letter.groups()
        wt = AA3_TO_AA1.get(wt3)
        mut = AA3_TO_AA1.get(mut3)
        if wt and mut:
            return MissenseVariant(position=int(pos), wt_aa=wt, mut_aa=mut)

    return None

