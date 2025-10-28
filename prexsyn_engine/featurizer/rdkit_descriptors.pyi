from collections.abc import Sequence

from .base import Featurizer

class RDKitDescriptorsFeaturizer(Featurizer):
    name: str
    num_evaluated_descriptors: int
    descriptor_names: Sequence[str]

    def __init__(
        self,
        name: str = "rdkit_descriptors",
        num_evaluated_descriptors: int = 4,
        descriptor_names: Sequence[str] = ...,
    ) -> None: ...
    def max_descriptor_index(self) -> int: ...
