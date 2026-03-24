import pytest


prexsyn_engine = pytest.importorskip("prexsyn_engine", exc_type=ImportError)
chemistry = prexsyn_engine.chemistry
Molecule = chemistry.Molecule
Reaction = chemistry.Reaction
Synthesis = chemistry.Synthesis
SynthesisError = chemistry.SynthesisError


REACTION_SMARTS = (
    "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]"
    ">>[NH:1]1[c:2][c:3][C:4](=[O:5])[N:7]([C:6])C1=O"
)
EXPECTED_PRODUCT_SMILES = "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O"


def make_test_reaction():
    return Reaction.from_smarts(REACTION_SMARTS, ["A", "B"])


def make_reactant_a():
    return Molecule.from_smiles("COC(=O)c1nscc1N")


def make_reactant_b():
    return Molecule.from_smiles("CN(C)CCCN")


def make_non_matching_reactant():
    return Molecule.from_smiles("CC")


def test_push_reaction_builds_expected_top_node_and_precursors():
    synthesis = Synthesis()
    reactant_a = make_reactant_a()
    reactant_b = make_reactant_b()
    reaction = make_test_reaction()

    synthesis.push_molecule(reactant_a)
    synthesis.push_molecule(reactant_b)
    synthesis.push_reaction(reaction, None)

    assert synthesis.stack_size() == 1
    top = synthesis.stack_top()
    assert top.size() == 1
    assert top.at(0).smiles() == EXPECTED_PRODUCT_SMILES

    precursors = top.precursors(0)
    assert len(precursors) == 2

    assert precursors[0].precursor_index == 0
    assert precursors[0].reactant_name == "B"
    assert precursors[0].item_index == 0
    assert precursors[0].molecule.smiles() == reactant_b.smiles()

    assert precursors[1].precursor_index == 1
    assert precursors[1].reactant_name == "A"
    assert precursors[1].item_index == 0
    assert precursors[1].molecule.smiles() == reactant_a.smiles()


def test_push_reaction_raises_when_stack_has_too_few_reactants():
    synthesis = Synthesis()
    synthesis.push_molecule(make_reactant_a())

    with pytest.raises(SynthesisError, match="Not enough reactants on the stack"):
        synthesis.push_reaction(make_test_reaction(), None)


def test_push_reaction_raises_when_no_products_are_produced():
    synthesis = Synthesis()
    synthesis.push_molecule(make_non_matching_reactant())
    synthesis.push_molecule(make_non_matching_reactant())

    with pytest.raises(SynthesisError, match="did not produce any products"):
        synthesis.push_reaction(make_test_reaction(), None)


def test_undo_restores_precursor_nodes_in_original_stack_order():
    synthesis = Synthesis()
    reactant_a = make_reactant_a()
    reactant_b = make_reactant_b()

    synthesis.push_molecule(reactant_a)
    synthesis.push_molecule(reactant_b)
    synthesis.push_reaction(make_test_reaction(), None)

    assert len(synthesis.nodes()) == 3
    assert synthesis.stack_size() == 1

    synthesis.undo()

    assert len(synthesis.nodes()) == 2
    assert synthesis.stack_size() == 2
    assert synthesis.stack_top(0).at(0).smiles() == reactant_b.smiles()
    assert synthesis.stack_top(1).at(0).smiles() == reactant_a.smiles()


def test_precursors_raises_on_invalid_item_index():
    synthesis = Synthesis()
    synthesis.push_molecule(make_reactant_a())
    synthesis.push_molecule(make_reactant_b())
    synthesis.push_reaction(make_test_reaction(), None)

    top = synthesis.stack_top()
    assert top.size() == 1

    with pytest.raises(IndexError):
        top.precursors(1)


def test_undo_raises_when_stack_is_empty():
    synthesis = Synthesis()

    with pytest.raises(SynthesisError, match="stack is empty"):
        synthesis.undo()
