import pathlib
import tempfile
from collections.abc import Iterable

import pytest
import rdkit.Chem

from prexsyn_engine import building_block_list


@pytest.fixture(scope="module")  # type: ignore[misc]
def temp_dir() -> Iterable[pathlib.Path]:
    with tempfile.TemporaryDirectory() as tmpdir:
        yield pathlib.Path(tmpdir)


def test_load_and_save(temp_dir: pathlib.Path) -> None:
    bb_list = building_block_list.BuildingBlockList.from_sdf("resources/test/chemspace_small_1/bb.sdf")
    bb_list.save(str(temp_dir / "test_building_blocks.bin"))

    bb_list2 = building_block_list.BuildingBlockList.load(str(temp_dir / "test_building_blocks.bin"))
    assert len(bb_list) == len(bb_list2)
    for i in range(min(len(bb_list), 10)):
        assert rdkit.Chem.MolToSmiles(bb_list[i]) == rdkit.Chem.MolToSmiles(bb_list2[i])
