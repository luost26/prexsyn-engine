import tempfile
from pathlib import Path

import pytest

from prexsyn_engine import chemspace, chemistry

Molecule = chemistry.Molecule
Reaction = chemistry.Reaction


def resource_path(name: str) -> Path:
    root = Path(__file__).resolve().parents[2]
    return root / "resources" / "test" / "chemspace_small_1" / name


def test_postfix_notation_basic_ops():
    pfn = chemspace.PostfixNotation()

    pfn.append(1, chemspace.PostfixNotationTokenType.BuildingBlock)

    t1 = chemspace.PostfixNotationToken()
    t1.index = 2
    t1.type = chemspace.PostfixNotationTokenType.BuildingBlock
    t2 = chemspace.PostfixNotationToken()
    t2.index = 3
    t2.type = chemspace.PostfixNotationTokenType.Reaction
    pfn.extend([t1, t2])

    tokens = pfn.tokens()
    assert len(tokens) == 3
    assert tokens[0].index == 1
    assert tokens[0].type == chemspace.PostfixNotationTokenType.BuildingBlock
    assert tokens[2].type == chemspace.PostfixNotationTokenType.Reaction

    pfn.pop_back()
    assert len(pfn.tokens()) == 2


def test_postfix_notation_invalid_type_raises():
    pfn = chemspace.PostfixNotation()

    with pytest.raises(TypeError):
        pfn.append(1, 99)


def test_building_block_library_add_get_and_serde_roundtrip():
    bb_lib = chemspace.BuildingBlockLibrary()

    entry = chemspace.BuildingBlockEntry()
    entry.molecule = Molecule.from_smiles("CCO")
    entry.identifier = "bb1"
    entry.labels = {"alcohol"}

    idx = bb_lib.add(entry)
    assert idx == 0
    assert bb_lib.size() == 1
    assert bb_lib.get(0).identifier == "bb1"
    assert bb_lib.get("bb1").molecule.smiles() == "CCO"

    with tempfile.NamedTemporaryFile() as tmp:
        bb_lib.serialize(tmp.name)
        cloned = chemspace.BuildingBlockLibrary.deserialize(tmp.name)
    assert cloned.size() == 1
    assert cloned.get(0).identifier == "bb1"


def test_building_block_library_duplicate_identifier_raises_specific_error():
    bb_lib = chemspace.BuildingBlockLibrary()

    entry = chemspace.BuildingBlockEntry()
    entry.molecule = Molecule.from_smiles("CC")
    entry.identifier = "dup"
    entry.labels = set()

    bb_lib.add(entry)
    with pytest.raises(
        chemspace.BuildingBlockLibraryError, match="duplicate identifier"
    ):
        bb_lib.add(entry)


def test_reaction_library_match_and_serde_roundtrip():
    rxn_lib = chemspace.ReactionLibrary()

    rxn_entry = chemspace.ReactionEntry()
    rxn_entry.reaction = Reaction.from_smarts(
        "[N:1].[O:2]>>[N:1][O:2]", ["amine", "alcohol"]
    )
    rxn_entry.name = "R1"
    rxn_lib.add(rxn_entry)

    matches = rxn_lib.match_reactants(Molecule.from_smiles("CCN"))
    assert len(matches) == 1
    assert matches[0].reaction_name == "R1"
    assert matches[0].reactant_name == "amine"

    with tempfile.NamedTemporaryFile() as tmp:
        rxn_lib.serialize(tmp.name)
        cloned = chemspace.ReactionLibrary.deserialize(tmp.name)
    assert cloned.size() == 1
    assert cloned.get(0).name == "R1"


def test_factory_loaders_from_test_resources():
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))

    assert bb_lib.size() > 0
    assert rxn_lib.size() > 0
    assert rxn_lib.get("ReactionA").name == "ReactionA"


def test_chemical_space_end_to_end_and_serde():
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))
    int_lib = chemspace.IntermediateLibrary()

    cs = chemspace.ChemicalSpace(bb_lib, rxn_lib, int_lib)
    cs.build_reactant_lists_for_building_blocks()
    assert cs.building_block_reactant_lists().num_matches() > 0

    cs.generate_intermediates()
    assert cs.int_lib().size() > 0

    cs.build_reactant_lists_for_intermediates()
    assert cs.intermediate_reactant_lists().num_matches() > 0

    rendered = cs.print_reactant_lists()
    assert "ReactionA" in rendered

    with tempfile.NamedTemporaryFile() as tmp:
        cs.serialize(tmp.name)
        stats = chemspace.ChemicalSpace.peek(tmp.name)
        assert stats.num_building_blocks == cs.bb_lib().size()
        assert stats.num_reactions == cs.rxn_lib().size()
        assert stats.num_intermediates == cs.int_lib().size()

        cloned = chemspace.ChemicalSpace.deserialize(tmp.name)
        assert cloned.bb_lib().size() == cs.bb_lib().size()
        assert cloned.rxn_lib().size() == cs.rxn_lib().size()
        assert cloned.int_lib().size() == cs.int_lib().size()


def test_chemspace_synthesis_add_and_undo():
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))
    int_lib = chemspace.IntermediateLibrary()
    cs = chemspace.ChemicalSpace(bb_lib, rxn_lib, int_lib)

    syn = cs.new_synthesis()

    r1 = syn.add_building_block("EN300-250786")
    r2 = syn.add_building_block("EN300-101318")
    r3 = syn.add_reaction("ReactionA", None)

    assert bool(r1) and r1.is_ok
    assert bool(r2) and r2.is_ok
    assert bool(r3) and r3.is_ok

    assert syn.count_building_blocks() == 2
    assert syn.count_reactions() == 1

    products = syn.products()
    assert len(products) > 0
    assert products[0].smiles() != ""

    undo_result = syn.undo()
    assert undo_result.is_ok


def test_chemspace_synthesis_add_building_block_accepts_int_and_string():
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))
    int_lib = chemspace.IntermediateLibrary()
    cs = chemspace.ChemicalSpace(bb_lib, rxn_lib, int_lib)

    syn = cs.new_synthesis()

    identifier = "EN300-250786"
    index = cs.bb_lib().get(identifier).index

    by_name = syn.add_building_block(identifier)
    by_index = syn.add_building_block(index)

    assert by_name.is_ok
    assert by_index.is_ok
    assert syn.count_building_blocks() == 2
