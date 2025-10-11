import pickle

import pytest
import rdkit.Chem
import rdkit.Chem.rdChemReactions

from prexsyn_engine import synthesis


@pytest.fixture(scope="function")  # type: ignore[misc]
def sample_1() -> synthesis.Synthesis:
    mol1 = rdkit.Chem.MolFromSmiles("COC(=O)c1nscc1N")
    mol2 = rdkit.Chem.MolFromSmiles("CN(C)CCCN")
    rxn = rdkit.Chem.rdChemReactions.ReactionFromSmarts(
        "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]>>" "[NH:1]1[c:2][c:3][C:4](=[O:5])[N:7]([C:6])C1=O"
    )
    rxn.Initialize()

    syn = synthesis.Synthesis()
    syn.push_mol(mol1)
    syn.push_mol(mol2)
    syn.push_reaction(rxn)
    return syn


def test_sample_1(sample_1: synthesis.Synthesis) -> None:
    products = sample_1.top().to_list()
    for product in products:
        print(rdkit.Chem.MolToSmiles(product))
    assert len(products) == 1
    assert rdkit.Chem.MolToSmiles(products[0]) == "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O"


def test_pfn_getitem(sample_1: synthesis.Synthesis) -> None:
    pfn = sample_1.get_postfix_notation()
    assert len(pfn) == 3
    assert isinstance(pfn[0], rdkit.Chem.Mol)
    assert isinstance(pfn[1], rdkit.Chem.Mol)
    assert isinstance(pfn[2], rdkit.Chem.rdChemReactions.ChemicalReaction)


def test_pickling_cpp(sample_1: synthesis.Synthesis) -> None:
    data = pickle.dumps(sample_1)
    new_syn: synthesis.Synthesis = pickle.loads(data)

    assert rdkit.Chem.MolToSmiles(new_syn.top().to_list()[0]) == "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O"
