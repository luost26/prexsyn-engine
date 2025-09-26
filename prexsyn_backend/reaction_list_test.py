import pathlib
import tempfile
from collections.abc import Iterable

import pytest
import rdkit.Chem.rdChemReactions

from . import reaction_list


@pytest.fixture(scope="module")  # type: ignore[misc]
def temp_dir() -> Iterable[pathlib.Path]:
    with tempfile.TemporaryDirectory() as tmpdir:
        yield pathlib.Path(tmpdir)


def test_load_and_save(temp_dir: pathlib.Path) -> None:
    rxn_list = reaction_list.ReactionList.from_txt("data/reactions/hartenfeller_button.txt")
    rxn_list.save(str(temp_dir / "test_reactions.bin"))

    rxn_list2 = reaction_list.ReactionList.load(str(temp_dir / "test_reactions.bin"))
    assert len(rxn_list) == len(rxn_list2)

    for i in range(len(rxn_list)):
        smarts1 = rdkit.Chem.rdChemReactions.ReactionToSmiles(rxn_list[i])
        smarts2 = rdkit.Chem.rdChemReactions.ReactionToSmiles(rxn_list2[i])
        assert smarts1 == smarts2
