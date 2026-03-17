import pickle

import pytest


prexsyn_engine = pytest.importorskip("prexsyn_engine", exc_type=ImportError)
chemistry = prexsyn_engine.chemistry
Molecule = chemistry.Molecule
Reaction = chemistry.Reaction
ReactionError = chemistry.ReactionError


REACTION_SMARTS = (
    "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]"
    ">>[NH:1]1[c:2][c:3][C:4](=[O:5])[N:7]([C:6])C1=O"
)
EXPECTED_PRODUCT_SMILES = "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O"


def make_reaction_with_named_reactants():
    return Reaction.from_smarts(REACTION_SMARTS, ["A", "B"])


def make_reactant_a():
    return Molecule.from_smiles("COC(=O)c1nscc1N")


def make_reactant_b():
    return Molecule.from_smiles("CN(C)CCCN")


def test_from_smarts_with_names_success():
    rxn = Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["amine", "alcohol"])

    assert rxn.num_reactants() == 2
    assert list(rxn.reactant_names()) == ["amine", "alcohol"]
    assert repr(rxn) == "<Reaction(amine, alcohol)>"


def test_from_smarts_auto_generates_reactant_names():
    rxn = Reaction.from_smarts("[C:1].[O:2]>>[C:1][O:2]")

    assert rxn.num_reactants() == 2
    assert list(rxn.reactant_names()) == ["R0", "R1"]
    assert repr(rxn) == "<Reaction(R0, R1)>"


def test_from_smarts_invalid_smarts_raises_reaction_error():
    with pytest.raises(ReactionError, match="Failed to parse SMARTS"):
        Reaction.from_smarts("not-a-smarts")


def test_from_smarts_mismatched_reactant_name_count_raises():
    with pytest.raises(
        ReactionError,
        match="Number of reactant names does not match number of reactant templates",
    ):
        Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["only_one_name"])


def test_from_smarts_duplicate_reactant_names_raises():
    with pytest.raises(ReactionError, match="Duplicate reactant name"):
        Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["dup", "dup"])


def test_match_reactants_returns_expected_matches():
    rxn = Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["amine", "alcohol"])
    mol = Molecule.from_smiles("CCN")

    matches = rxn.match_reactants(mol)

    assert isinstance(matches, list)
    assert len(matches) == 1
    assert matches[0]["index"] == 0
    assert matches[0]["name"] == "amine"
    assert matches[0]["count"] >= 1


def test_match_reactants_none_molecule_raises_reaction_error():
    rxn = Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["amine", "alcohol"])

    with pytest.raises(ReactionError, match="Molecule pointer is null"):
        rxn.match_reactants(None)


def test_pickle_roundtrip_preserves_reaction_metadata():
    rxn = Reaction.from_smarts("[N:1].[O:2]>>[N:1][O:2]", ["amine", "alcohol"])

    dumped = pickle.dumps(rxn)
    loaded = pickle.loads(dumped)

    assert loaded.num_reactants() == 2
    assert list(loaded.reactant_names()) == ["amine", "alcohol"]
    assert repr(loaded) == "<Reaction(amine, alcohol)>"


def test_apply_dict_success_returns_expected_product():
    rxn = make_reaction_with_named_reactants()

    outcomes = rxn.apply_dict({"A": make_reactant_a(), "B": make_reactant_b()})

    assert len(outcomes) == 1
    assert outcomes[0].empty() is False
    assert outcomes[0].num_products() == 1
    assert outcomes[0].main_product().smiles() == EXPECTED_PRODUCT_SMILES


def test_apply_dict_unknown_reactant_name_raises_reaction_error():
    rxn = make_reaction_with_named_reactants()

    with pytest.raises(ReactionError, match="Unknown reactant name"):
        rxn.apply_dict({"A": make_reactant_a(), "C": make_reactant_b()})


def test_apply_list_success_tracks_assignment_and_product():
    rxn = make_reaction_with_named_reactants()

    outcomes = rxn.apply_list([make_reactant_a(), make_reactant_b()])

    assert len(outcomes) == 1
    assert outcomes[0].reactant_names == ["A", "B"]
    assert outcomes[0].num_products() == 1
    assert outcomes[0].main_product().smiles() == EXPECTED_PRODUCT_SMILES


def test_apply_list_wrong_reactant_count_raises_reaction_error():
    rxn = make_reaction_with_named_reactants()

    with pytest.raises(
        ReactionError,
        match="Number of reactants provided does not match number of reactant templates",
    ):
        rxn.apply_list([make_reactant_a()])


def test_apply_dict_can_return_multiple_outcomes():
    # Each ethane reactant has two matching carbon atoms, yielding multiple outcomes.
    rxn = Reaction.from_smarts("[C:1].[C:2]>>[C:1][C:2]", ["A", "B"])

    outcomes = rxn.apply_dict(
        {"A": Molecule.from_smiles("CC"), "B": Molecule.from_smiles("CC")}
    )

    assert len(outcomes) > 1
    assert all(outcome.num_products() >= 1 for outcome in outcomes)


def test_apply_list_can_return_multiple_outcomes_from_assignments():
    # Symmetric templates with identical reactants produce outcomes for multiple assignments.
    rxn = Reaction.from_smarts("[C:1].[C:2]>>[C:1][C:2]", ["A", "B"])

    outcomes = rxn.apply_list([Molecule.from_smiles("C"), Molecule.from_smiles("C")])

    assert len(outcomes) > 1
    assert ["A", "B"] in [outcome.reactant_names for outcome in outcomes]
    assert ["B", "A"] in [outcome.reactant_names for outcome in outcomes]
