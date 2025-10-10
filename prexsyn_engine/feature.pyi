import abc
from typing import Any, TypeAlias

import numpy as np

FeatureDict: TypeAlias = dict[str, int | float | str | bool | np.ndarray[Any, Any]]

class FeatureBuilder(abc.ABC): ...

class PyDictBuilder:
    def __init__(self) -> None: ...
    def get(self) -> FeatureDict: ...
