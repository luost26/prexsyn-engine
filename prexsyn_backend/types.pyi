from collections.abc import MutableSequence
from typing import TypeAlias

from rdkit.Chem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction as Reaction

class Path:
    def __init__(self, _: str, /) -> None: ...

MolVector: TypeAlias = MutableSequence[Mol]
ReactionVector: TypeAlias = MutableSequence[Reaction]

class MolSet:
    def to_list(self) -> list[Mol]: ...
    def __len__(self) -> int: ...

__all__ = ["Mol", "Reaction", "Path", "MolVector", "ReactionVector", "MolSet"]
