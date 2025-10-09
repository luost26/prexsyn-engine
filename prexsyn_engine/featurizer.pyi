import abc
from collections.abc import Sequence
from typing import Any, Self, TypeAlias

import numpy as np

from .synthesis import Synthesis

FeatureDict: TypeAlias = dict[str, int | float | str | bool | np.ndarray[Any, Any]]

class FeatureBuilder(abc.ABC): ...

class PyDictBuilder:
    def __init__(self) -> None: ...
    def get(self) -> FeatureDict: ...
    def erase_type(self) -> FeatureBuilder: ...

class Featurizer(abc.ABC):
    @abc.abstractmethod
    def __call__(self, synthesis: Synthesis, builder: FeatureBuilder) -> None: ...

class FeaturizerSet(Featurizer):
    def __init__(self) -> None: ...
    def __call__(self, synthesis: Synthesis, builder: FeatureBuilder) -> None: ...
    def add(self, featurizer: Featurizer) -> Self: ...

class ProductStructureFeaturizerOption:
    fp_types: list[str]
    embedder_name_template: str
    scaffold: bool
    scaffold_fp_type: str
    scaffold_embedder_name: str
    num_fragments: int
    fragment_fp_type: str
    fragment_embedder_name: str

class ProductStructureFeaturizer(Featurizer):
    def __init__(
        self,
        option: ProductStructureFeaturizerOption = ProductStructureFeaturizerOption(),
    ) -> None: ...
    def __call__(self, synthesis: Synthesis, builder: FeatureBuilder) -> None: ...

RDKitPropertyList: TypeAlias = Sequence[str]

class ProductRDKitPropertyFeaturizerOption:
    name: str
    num_evaluated_properties: int
    rdkit_property_index_offset: int
    rdkit_properties: RDKitPropertyList

class ProductRDKitPropertyFeaturizer(Featurizer):
    def __init__(
        self,
        option: ProductRDKitPropertyFeaturizerOption = ProductRDKitPropertyFeaturizerOption(),
    ) -> None: ...
    def __call__(self, synthesis: Synthesis, builder: FeatureBuilder) -> None: ...
    def max_property_index(self) -> int: ...

class PostfixNotationTokenDef:
    PAD: int
    END: int
    START: int
    BB: int
    RXN: int

class PostfixNotationFeaturizerOption:
    length: int
    token_def: PostfixNotationTokenDef

class PostfixNotationFeaturizer(Featurizer):
    option: PostfixNotationFeaturizerOption
    def __init__(
        self,
        option: PostfixNotationFeaturizerOption = PostfixNotationFeaturizerOption(),
    ) -> None: ...
    def __call__(self, synthesis: Synthesis, builder: FeatureBuilder) -> None: ...
