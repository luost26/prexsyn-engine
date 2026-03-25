import pickle

import pytest
import prexsyn_engine


Chem = pytest.importorskip("rdkit.Chem", exc_type=ImportError)

chemistry = prexsyn_engine.chemistry
Molecule = chemistry.Molecule
MoleculeError = chemistry.MoleculeError


def test_from_smiles_success():
    mol = Molecule.from_smiles("CCO")

    assert mol.smiles() == "CCO"
    assert mol.num_heavy_atoms() == 3


def test_from_smiles_invalid_raises_molecule_error():
    with pytest.raises(MoleculeError, match="Failed to parse SMILES"):
        Molecule.from_smiles("not-a-smiles")


def test_repr_contains_smiles():
    mol = Molecule.from_smiles("CO")

    assert repr(mol) == "<Molecule CO>"


def test_pickle_roundtrip_preserves_smiles_and_properties():
    mol = Molecule.from_smiles("CCN")

    dumped = pickle.dumps(mol)
    loaded = pickle.loads(dumped)

    assert loaded.smiles() == "CCN"
    assert loaded.num_heavy_atoms() == 3


def test_largest_fragment_returns_main_component():
    # Ethane + water separated by dot notation; ethane is the largest fragment.
    mol = Molecule.from_smiles("CC.O")

    largest = mol.largest_fragment()

    assert largest.smiles() == "CC"
    assert largest.num_heavy_atoms() == 2


def test_from_rdkit_mol_creates_prexsyn_molecule():
    rdk_mol = Chem.MolFromSmiles("CCO")

    mol = Molecule.from_rdkit_mol(rdk_mol)

    assert mol.smiles() == "CCO"
    assert mol.num_heavy_atoms() == 3


def test_to_rdkit_mol_returns_rdkit_molecule():
    mol = Molecule.from_smiles("CCN")

    rdk_mol = mol.to_rdkit_mol()

    assert isinstance(rdk_mol, Chem.rdchem.Mol)
    assert Chem.MolToSmiles(rdk_mol) == "CCN"


def test_rdkit_conversion_roundtrip_preserves_canonical_smiles():
    input_smiles = "OC1=CC=CC=C1"
    canonical_input = Chem.MolToSmiles(Chem.MolFromSmiles(input_smiles))

    mol = Molecule.from_rdkit_mol(Chem.MolFromSmiles(input_smiles))
    rdk_mol = mol.to_rdkit_mol()
    canonical_output = Chem.MolToSmiles(rdk_mol)

    assert canonical_output == canonical_input
