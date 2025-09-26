from collections.abc import Sequence

from .synthesis import Synthesis
from .types import Mol, Reaction

def get_suitable_reactant_indices_for_mol(reaction: Reaction, mol: Mol) -> Sequence[int]: ...
def get_suitable_reactant_indices_for_synthesis(reaction: Reaction, synthesis: Synthesis) -> Sequence[int]: ...
