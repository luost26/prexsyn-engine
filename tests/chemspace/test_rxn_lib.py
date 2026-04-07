import json
from pathlib import Path

import pytest

from prexsyn_engine import chemistry, chemspace


def resource_path(name: str) -> Path:
    root = Path(__file__).resolve().parents[2]
    return root / "resources" / "test" / name


def test_rxn_lib_from_json_loads_example_file():
    rxn_lib = chemspace.rxn_lib_from_json(resource_path("reaction.json"))

    assert rxn_lib.size() == 3

    suzuki = rxn_lib.get("Suzuki Coupling")
    assert suzuki.name == "Suzuki Coupling"
    assert suzuki.reaction.num_reactants() == 2
    assert set(suzuki.reaction.reactant_names()) == {"Halides", "Boronates"}

    amide = rxn_lib.get("Amide Coupling")
    assert amide.reaction.num_reactants() == 2
    assert set(amide.reaction.reactant_names()) == {"Acids", "Amines"}

    alkynylation = rxn_lib.get("Alkynylation")
    assert alkynylation.reaction.num_reactants() == 1
    assert alkynylation.reaction.reactant_names() == ["Halides"]


def test_rxn_lib_from_json_reactions_apply_to_concrete_molecules():
    rxn_lib = chemspace.rxn_lib_from_json(resource_path("reaction.json"))

    suzuki = rxn_lib.get("Suzuki Coupling").reaction
    suzuki_outcomes = suzuki.apply(
        {
            "Halides": chemistry.Molecule.from_smiles("c1ccccc1Br"),
            "Boronates": chemistry.Molecule.from_smiles("OB(O)c1ccccc1"),
        }
    )
    assert len(suzuki_outcomes) > 0
    assert suzuki_outcomes[0].num_products() > 0
    assert suzuki_outcomes[0].main_product().smiles() != ""

    amide = rxn_lib.get("Amide Coupling").reaction
    amide_outcomes = amide.apply(
        {
            "Acids": chemistry.Molecule.from_smiles("CC(=O)O"),
            "Amines": chemistry.Molecule.from_smiles("CN"),
        }
    )
    assert len(amide_outcomes) > 0
    assert amide_outcomes[0].num_products() > 0
    assert "N" in amide_outcomes[0].main_product().smiles()
    assert "=O" in amide_outcomes[0].main_product().smiles()

    alkynylation = rxn_lib.get("Alkynylation").reaction
    alkynylation_outcomes = alkynylation.apply(
        {
            "Halides": chemistry.Molecule.from_smiles("c1ccccc1Br"),
        }
    )
    assert len(alkynylation_outcomes) > 0
    assert alkynylation_outcomes[0].num_products() > 0
    assert "#" in alkynylation_outcomes[0].main_product().smiles()


def test_rxn_lib_from_json_raises_when_ignore_errors_false(tmp_path: Path):
    path = tmp_path / "invalid_reactions.json"
    path.write_text(
        json.dumps(
            [
                {
                    "name": "Broken",
                    "reactants": {"A": "[C:1]-[Br]"},
                    # missing required "product" field
                }
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(Exception):
        chemspace.rxn_lib_from_json(path, ignore_errors=False)


def test_rxn_lib_from_json_skips_invalid_entries_when_ignore_errors_true(
    tmp_path: Path,
):
    path = tmp_path / "mixed_reactions.json"
    path.write_text(
        json.dumps(
            [
                {
                    "name": "Broken",
                    "reactants": {"A": "[C:1]-[Br]"},
                    # missing required "product" field
                },
                {
                    "name": "Valid",
                    "reactants": {"A": "[C:1]-[Br]"},
                    "product": "[C:1]-[C]",
                },
            ]
        ),
        encoding="utf-8",
    )

    rxn_lib = chemspace.rxn_lib_from_json(path, ignore_errors=True)

    assert rxn_lib.size() == 1
    assert rxn_lib.get("Valid").name == "Valid"
    with pytest.raises(Exception):
        rxn_lib.get("Broken")
