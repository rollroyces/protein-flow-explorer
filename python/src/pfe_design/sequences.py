"""Residue-type indices and sequences.

Two things need to agree: the integer residue types the network consumes, and
the single-letter sequence a structure predictor such as ESMFold wants.

The convention is deliberately small:

* index ``0`` is "unknown" and is written as ``X``;
* indices ``1..20`` are the twenty standard amino acids, in the order below.

Keeping the two in one place is what lets the packing backends receive a
sequence rather than a list of integers.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

__all__ = [
    "AA1_ORDER",
    "AA1_TO_INDEX",
    "INDEX_TO_AA1",
    "AA3_TO_AA1",
    "sequence_from_residue_types",
    "residue_types_from_sequence",
    "canonical_sequence",
]

# The twenty standard amino acids, in the order the residue-type indices use.
AA1_ORDER = "ARNDCQEGHILKMFPSTWYV"

AA1_TO_INDEX = {aa: index + 1 for index, aa in enumerate(AA1_ORDER)}
INDEX_TO_AA1 = {index + 1: aa for index, aa in enumerate(AA1_ORDER)}

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

# Modified residues that a structure file may carry, mapped onto the residue
# they are derived from.  A predictor cannot be asked about MSE; it can be
# asked about MET.
MODIFIED_AA3 = {
    "MSE": "M",  # selenomethionine
    "SEC": "C",  # selenocysteine
    "CSU": "C",  # cysteine, alternate
    "CSX": "C",  # cysteine, unknown state
    "PYL": "K",  # pyrrolysine
    "HYP": "P",  # hydroxyproline
    "MLY": "K",  # methyllysine
    "SEP": "S",  # phosphoserine
    "TPO": "T",  # phosphothreonine
    "PTR": "Y",  # phosphotyrosine
    "HIC": "H",  # histidine, alternate
    "PCA": "E",  # pyroglutamate
    "CGU": "E",  # glutamate, alternate
    "LLP": "K",  # lysine, alternate
}

# Ambiguous or unknown residues have no single-letter equivalent.
AMBIGUOUS_AA3 = {
    "ASX": "X",
    "GLX": "X",
    "UNK": "X",
    "XAA": "X",
    "UNL": "X",
}


def sequence_from_residue_types(residue_types: Iterable[int]) -> str:
    """Write a single-letter sequence from residue-type indices.

    An index of ``0``, or one outside the known range, becomes ``X``: an
    unknown residue is written as unknown rather than guessed at.
    """
    letters = []
    for value in residue_types:
        index = int(value)
        letters.append(INDEX_TO_AA1.get(index, "X"))
    return "".join(letters)


def residue_types_from_sequence(sequence: str, format: str = "auto") -> list[int]:
    """Read residue-type indices from a sequence.

    Three-letter codes are accepted, which is what a structure file gives,
    and the modified residues in :data:`MODIFIED_AA3` are mapped onto the
    residue they are derived from.

    The ``format`` argument says how to read the string:

    ``"auto"`` (the default)
        A string with more than one whitespace-separated part is read as
        three-letter codes; a single part is read as one letter per
        residue.  The rule is not cosmetic: every three-letter code is also
        a valid run of three one-letter codes, so ``"MET"`` is ambiguous,
        and reading it one way would silently change the length of the
        sequence - which is the sequence a structure predictor is given.
    ``"one-letter"``
        One residue per character, always.
    ``"three-letter"``
        One residue per whitespace-separated part; a part that is not a known
        code is read as one letter per character.
    """
    indices: list[int] = []
    for token in _tokenise(sequence.strip(), format):
        letter = _to_single_letter(token)
        indices.append(AA1_TO_INDEX.get(letter, 0))
    return indices


def canonical_sequence(sequence: str, format: str = "auto") -> str:
    """Normalise a sequence to single-letter codes.

    Modified residues become the residue they are derived from and ambiguous
    ones become ``X``.
    """
    return sequence_from_residue_types(residue_types_from_sequence(sequence, format))


def _tokenise(text: str, format: str = "auto") -> list[str]:
    """Split a sequence into residue tokens.

    The hard part is that every three-letter code is also a valid run of
    three one-letter codes, so a string like ``MET`` is ambiguous: read as
    three-letter codes it is one residue, read as one-letter codes it is
    three.  Guessing wrong changes the length of the sequence, which is the
    sequence a predictor is given, so the guess is made explicitly:

    * more than one whitespace-separated part means three-letter codes;
    * a single part means one letter per residue.

    ``format`` overrides both.
    """
    parts = [part for part in re.split(r"[\s,]+", text) if part]
    if not parts:
        return []

    if format == "one-letter":
        return list("".join(parts).upper())
    if format == "three-letter":
        tokens: list[str] = []
        for part in parts:
            tokens.extend(_three_letter_tokens(part))
        return tokens
    if format != "auto":
        raise ValueError(f"unknown sequence format: {format!r}")

    if len(parts) > 1:
        tokens = []
        for part in parts:
            tokens.extend(_three_letter_tokens(part))
        return tokens
    return list("".join(parts).upper())


def _three_letter_tokens(part: str) -> list[str]:
    """Split one whitespace-separated part into residue tokens."""
    upper = part.upper()
    known = {*AA3_TO_AA1, *MODIFIED_AA3, *AMBIGUOUS_AA3}
    if len(upper) == 3 and upper in known:
        return [upper]
    if len(upper) % 3 == 0 and len(upper) > 0:
        chunks = [upper[i : i + 3] for i in range(0, len(upper), 3)]
        if all(chunk in known for chunk in chunks):
            return chunks
    return list(upper)


def _to_single_letter(token: str) -> str:
    if len(token) == 3:
        if token in AMBIGUOUS_AA3:
            return AMBIGUOUS_AA3[token]
        if token in MODIFIED_AA3:
            return MODIFIED_AA3[token]
        return AA3_TO_AA1.get(token, "X")
    return token if token in AA1_TO_INDEX else "X"
