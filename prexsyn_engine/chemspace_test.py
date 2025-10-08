import pathlib
import tempfile
from collections.abc import Iterable

import pytest
import rdkit.Chem

from . import chemspace


@pytest.fixture  # type: ignore[misc]
def temp_dir() -> Iterable[pathlib.Path]:
    with tempfile.TemporaryDirectory() as tmpdir:
        yield pathlib.Path(tmpdir)


def test_create_chemspace(temp_dir: pathlib.Path) -> None:
    # Test the creation of a chemspace
    csd = (
        chemspace.ChemicalSpaceDefinitionBuilder()
        .building_blocks_from_sdf("data/building_blocks/mcule_subset.sdf")
        .reactions_from_txt("data/reactions/hartenfeller_button.txt")
        .secondary_building_blocks_from_single_reaction()
        .build_primary_index()
        .build_secondary_index()
        .build()
    )

    csd.save(str(temp_dir / "test_chemspace.bin"))
    chemspace.ChemicalSpaceDefinitionBuilder().all_from_cache(str(temp_dir / "test_chemspace.bin")).build()


def test_generate_synthesis() -> None:
    csd = (
        chemspace.ChemicalSpaceDefinitionBuilder()
        .building_blocks_from_sdf("data/building_blocks/mcule_subset.sdf")
        .reactions_from_txt("data/reactions/hartenfeller_button.txt")
        .no_secondary_building_blocks()
        .build_primary_index()
        .build_secondary_index()
        .build()
    )
    generator = chemspace.SynthesisGenerator(csd)
    for _ in range(10):
        syn = generator.next()
        print(
            syn.count_building_blocks(),
            syn.count_reactions(),
            [rdkit.Chem.MolToSmiles(m) for m in syn.top().to_list()],
        )
