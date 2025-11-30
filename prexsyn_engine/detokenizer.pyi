from typing import Any

import numpy as np

from . import building_block_list, featurizer, reaction_list, synthesis

class Detokenizer:
    def __init__(
        self,
        building_blocks: building_block_list.BuildingBlockList,
        reactions: reaction_list.ReactionList,
        token_def: featurizer.synthesis.PostfixNotationTokenDef = ...,
    ) -> None: ...
    def __call__(
        self,
        token_types: np.ndarray[Any, Any],
        bb_indices: np.ndarray[Any, Any],
        rxn_indices: np.ndarray[Any, Any],
    ) -> synthesis.SynthesisVector: ...
