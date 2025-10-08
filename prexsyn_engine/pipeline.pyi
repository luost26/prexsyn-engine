from typing import Any

import numpy as np

from .chemspace import ChemicalSpaceDefinition, SynthesisGeneratorOption
from .featurizer import FeaturizerSet, PyDictBuilder
from .synthesis import Synthesis

class DataPipeline:
    def __init__(
        self,
        num_threads: int,
        csd: ChemicalSpaceDefinition,
        gen_option: SynthesisGeneratorOption,
        featurizer: FeaturizerSet,
    ) -> None: ...
    def start(self) -> None: ...
    def get(self) -> tuple[Synthesis, PyDictBuilder]: ...
    def stop(self) -> None: ...

class DataPipelineV2:
    def __init__(
        self,
        num_threads: int,
        csd: ChemicalSpaceDefinition,
        gen_option: SynthesisGeneratorOption,
        featurizer: FeaturizerSet,
        base_seed: int = ...,
    ) -> None: ...
    def start(self) -> None: ...
    def get(self, n: int) -> dict[str, np.ndarray[Any, Any]]: ...
    def stop(self) -> None: ...
