from collections.abc import MutableMapping
from typing import Any, TypeAlias

import numpy as np
from rdkit.Chem import BondType
from rdkit.Chem.rdMolChemicalFeatures import (
    MolChemicalFeature,
    MolChemicalFeatureFactory,
)

from .types import Mol

class PharmacophoreNode:
    family: str
    type: str

class PharmacophoreGraph:
    def set_edge(self, u: int, v: int, weight: float = ...) -> None: ...
    def unset_edge(self, u: int, v: int) -> None: ...
    def has_edge(self, u: int, v: int) -> bool: ...
    def get_nodes(self) -> list[PharmacophoreNode]: ...
    def get_edges(self) -> list[tuple[int, int, float]]: ...
    def get_min_spanning_tree(self) -> "PharmacophoreGraph": ...
    def get_subgraph_with_n_edges(self, n_edges: int) -> "PharmacophoreGraph": ...

BondWeights: TypeAlias = MutableMapping[BondType, float]

bond_weights_EMPTY: BondWeights
bond_weights_RELATIVE_BOND_LENGTH: BondWeights

def create_pharmacophore_graph(
    mol: Mol,
    features: Any | MolChemicalFeatureFactory | list[MolChemicalFeature],
    bond_weights: BondWeights = bond_weights_EMPTY,
    default_bond_weight: float = 1.0,
) -> PharmacophoreGraph: ...
def subgraph_similarity(
    g1: PharmacophoreGraph,
    g2: PharmacophoreGraph,
    subgraph_nodes_1: list[int],
    subgraph_nodes_2: list[int],
    lambda_weight: float = 1.0,
) -> tuple[float, float]: ...
def pharmacophore_similarity(
    g1: PharmacophoreGraph,
    g2: PharmacophoreGraph,
    lambda_weight: float = 1.0,
) -> tuple[float, float]: ...
def pairwise_pharmacophore_similarity(
    graphs1: list[PharmacophoreGraph],
    graphs2: list[PharmacophoreGraph],
    lambda_weight: float = 1.0,
) -> tuple[np.ndarray[tuple[int, int], np.dtype[np.float32]], np.ndarray[tuple[int, int], np.dtype[np.float32]]]: ...
