import prexsyn_engine.chemistry
import prexsyn_engine.chemspace
import typing

class EnumeratorConfig:
    heavy_atom_limit: int
    max_building_blocks: int
    max_outcomes_per_reaction: int
    selectivity_cutoff: int
    def __init__(self) -> None: ...

class RandomEnumerator:
    def __init__(self, chemical_space: prexsyn_engine.chemspace.ChemicalSpace, config: EnumeratorConfig = ..., random_seed: typing.SupportsInt | typing.SupportsIndex | None = ...) -> None: ...
    def next(self) -> prexsyn_engine.chemspace.Synthesis: ...
    def next_with_product(self) -> tuple[prexsyn_engine.chemspace.Synthesis, prexsyn_engine.chemistry.Molecule]: ...
