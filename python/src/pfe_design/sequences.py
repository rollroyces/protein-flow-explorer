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


def residue_types_from_sequence(sequence: str) -> list[int]:
    """Read residue-type indices from a single-letter sequence.

    Three-letter codes are accepted too, which is what a structure file
    gives, and the modified residues in :data:`MODIFIED_AA3` are mapped onto
    the residue they are derived from.
    """
    indices: list[int] = []
    text = sequence.strip()
    for token in _tokenise(text):
        letter = _to_single_letter(token)
        indices.append(AA1_TO_INDEX.get(letter, 0))
    return indices


def canonical_sequence(sequence: str) -> str:
    """Normalise a sequence to single-letter codes.

    Modified residues become the residue they are derived from and ambiguous
    ones become ``X``.
    """
    return sequence_from_residue_types(residue_types_from_sequence(sequence))


def _tokenise(text: str) -> list[str]:
    """Split a sequence into residue tokens.

    A three-letter code is one token; a run of single-letter codes is one
    token per letter.  Whitespace is ignored, so a sequence that was written
    out in columns still reads correctly.
    """
    tokens: list[str] = []
    index = 0
    upper = text.upper()
    while index < len(upper):
        character = upper[index]
        if character.isspace():
            index += 1
            continue
        chunk = upper[index : index + 3]
        if len(chunk) == 3 and chunk in {*AA3_TO_AA1, *MODIFIED_AA3, *AMBIGUOUS_AA3}:
            tokens.append(chunk)
            index += 3
            continue
        tokens.append(character)
        index += 1
    return tokens


def _to_single_letter(token: str) -> str:
    if len(token) == 3:
        if token in AMBIGUOUS_AA3:
            return AMBIGUOUS_AA3[token]
        if token in MODIFIED_AA3:
            return MODIFIED_AA3[token]
        return AA3_TO_AA1.get(token, "X")
    return token if token in AA1_TO_INDEX else "X"
