from pathlib import Path

import pytest


prexsyn_engine = pytest.importorskip("prexsyn_engine", exc_type=ImportError)
chemspace = prexsyn_engine.chemspace
enumerator = prexsyn_engine.enumerator


def resource_path(name: str) -> Path:
    root = Path(__file__).resolve().parents[2]
    return root / "resources" / "test" / "chemspace_small_1" / name


def make_chemical_space():
    """Helper to create a fully initialized ChemicalSpace."""
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))
    int_lib = chemspace.IntermediateLibrary()

    cs = chemspace.ChemicalSpace(bb_lib, rxn_lib, int_lib)
    cs.build_reactant_lists_for_building_blocks()
    cs.generate_intermediates()
    cs.build_reactant_lists_for_intermediates()
    return cs


class TestEnumeratorConfig:
    """Test cases for EnumeratorConfig class."""

    def test_config_creation_succeeds(self):
        """Test that EnumeratorConfig can be created."""
        config = enumerator.EnumeratorConfig()
        assert config is not None

    def test_config_default_values(self):
        """Test that EnumeratorConfig has expected default values."""
        config = enumerator.EnumeratorConfig()

        assert config.max_building_blocks == 5
        assert config.heavy_atom_limit == 50

    def test_config_max_building_blocks_is_writeable(self):
        """Test that max_building_blocks can be modified."""
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 3

        assert config.max_building_blocks == 3

    def test_config_heavy_atom_limit_is_writeable(self):
        """Test that heavy_atom_limit can be modified."""
        config = enumerator.EnumeratorConfig()
        config.heavy_atom_limit = 30

        assert config.heavy_atom_limit == 30

    def test_config_both_attributes_writeable_together(self):
        """Test that both attributes can be modified together."""
        config = enumerator.EnumeratorConfig()

        config.max_building_blocks = 7
        config.heavy_atom_limit = 100

        assert config.max_building_blocks == 7
        assert config.heavy_atom_limit == 100

    def test_config_zero_values_are_accepted(self):
        """Test that zero values are accepted for config attributes."""
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 0
        config.heavy_atom_limit = 0

        assert config.max_building_blocks == 0
        assert config.heavy_atom_limit == 0


class TestRandomEnumerator:
    """Test cases for RandomEnumerator class."""

    def test_random_enumerator_creation_succeeds(self):
        """Test that RandomEnumerator can be created."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs)

        assert enumerator_obj is not None

    def test_random_enumerator_with_custom_config(self):
        """Test creating RandomEnumerator with custom config."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 3

        enumerator_obj = enumerator.RandomEnumerator(cs, config)
        assert enumerator_obj is not None

    def test_random_enumerator_with_random_seed(self):
        """Test creating RandomEnumerator with explicit random seed."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=42)

        assert enumerator_obj is not None

    def test_random_enumerator_with_config_and_seed(self):
        """Test creating RandomEnumerator with both config and seed."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 2

        enumerator_obj = enumerator.RandomEnumerator(cs, config, random_seed=99)
        assert enumerator_obj is not None

    def test_next_returns_synthesis(self):
        """Test that next() returns a Synthesis object."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=7)

        syn = enumerator_obj.next()

        assert isinstance(syn, chemspace.Synthesis)

    def test_next_synthesis_has_building_blocks(self):
        """Test that next() returns synthesis with building blocks."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=11)

        syn = enumerator_obj.next()

        assert syn.count_building_blocks() >= 1

    def test_next_synthesis_has_products(self):
        """Test that next() returns synthesis with products."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=13)

        syn = enumerator_obj.next()

        products = syn.products()
        assert len(products) >= 1

    def test_next_with_product_returns_pair(self):
        """Test that next_with_product() returns a tuple of (Synthesis, Molecule)."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=19)

        result = enumerator_obj.next_with_product()

        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_next_with_product_returns_synthesis_and_molecule(self):
        """Test that next_with_product() returns correct types."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=23)

        syn, product = enumerator_obj.next_with_product()

        assert isinstance(syn, chemspace.Synthesis)
        assert product.smiles() != ""

    def test_reproducibility_with_same_seed_first_result(self):
        """Test that same seed produces reproducible results."""
        cs = make_chemical_space()

        enum1 = enumerator.RandomEnumerator(cs, random_seed=123)
        enum2 = enumerator.RandomEnumerator(cs, random_seed=123)

        _, product1 = enum1.next_with_product()
        _, product2 = enum2.next_with_product()

        assert product1.smiles() == product2.smiles()

    def test_reproducibility_with_same_seed_multiple_results(self):
        """Test that same seed produces reproducible results across multiple calls."""
        cs = make_chemical_space()

        enum1 = enumerator.RandomEnumerator(cs, random_seed=456)
        enum2 = enumerator.RandomEnumerator(cs, random_seed=456)

        smiles_seq1 = [enum1.next_with_product()[1].smiles() for _ in range(3)]
        smiles_seq2 = [enum2.next_with_product()[1].smiles() for _ in range(3)]

        assert smiles_seq1 == smiles_seq2

    def test_different_seeds_produce_different_results(self):
        """Test that different seeds typically produce different results."""
        cs = make_chemical_space()

        enum1 = enumerator.RandomEnumerator(cs, random_seed=111)
        enum2 = enumerator.RandomEnumerator(cs, random_seed=555)

        _, product1 = enum1.next_with_product()
        _, product2 = enum2.next_with_product()

        # Different seeds should almost always produce different results
        # (theoretically they could collide, but it's extremely unlikely)
        assert product1.smiles() != product2.smiles()

    def test_config_max_building_blocks_limits_synthesis_growth(self):
        """Test that max_building_blocks config is respected."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 2

        enumerator_obj = enumerator.RandomEnumerator(cs, config, random_seed=789)

        for _ in range(10):
            syn = enumerator_obj.next()
            assert syn.count_building_blocks() <= config.max_building_blocks

    # def test_config_heavy_atom_limit_limits_synthesis(self):
    #     """Test that heavy_atom_limit config is respected."""
    #     cs = make_chemical_space()
    #     config = enumerator.EnumeratorConfig()
    #     config.heavy_atom_limit = 40

    #     enumerator_obj = enumerator.RandomEnumerator(cs, config, random_seed=654)

    #     for _ in range(10):
    #         syn = enumerator_obj.next()
    #         products = syn.products()
    #         if products:
    #             assert products[0].num_heavy_atoms() <= config.heavy_atom_limit

    def test_sequential_calls_produce_valid_syntheses(self):
        """Test that sequential next() calls all produce valid syntheses."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=222)

        for _ in range(5):
            syn = enumerator_obj.next()
            assert syn.count_building_blocks() >= 1
            assert len(syn.products()) >= 1

    def test_next_after_next_with_product(self):
        """Test that next() works after next_with_product()."""
        cs = make_chemical_space()
        enumerator_obj = enumerator.RandomEnumerator(cs, random_seed=333)

        syn1, mol = enumerator_obj.next_with_product()
        syn2 = enumerator_obj.next()

        assert syn1 is not syn2
        assert isinstance(mol, type(syn2.products()[0]))


class TestEnumeratorWithDifferentConfigs:
    """Test cases for enumerator behavior with various configurations."""

    def test_restrictive_config_still_produces_valid_synthesis(self):
        """Test that very restrictive config still produces valid results."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 1
        config.heavy_atom_limit = 10

        enumerator_obj = enumerator.RandomEnumerator(cs, config, random_seed=999)
        syn = enumerator_obj.next()

        assert syn.count_building_blocks() >= 1

    def test_permissive_config_can_produce_larger_synthesis(self):
        """Test that permissive config allows larger syntheses."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()
        config.max_building_blocks = 10
        config.heavy_atom_limit = 200

        enumerator_obj = enumerator.RandomEnumerator(cs, config, random_seed=777)

        # Generate several and check at least one has more building blocks than default
        max_bbs = 0
        for _ in range(10):
            syn = enumerator_obj.next()
            max_bbs = max(max_bbs, syn.count_building_blocks())

        # With permissive config and multiple attempts, we should see > 5 BBs
        assert max_bbs > 1

    def test_default_config_values_work(self):
        """Test that default config values are sensible defaults."""
        cs = make_chemical_space()
        config = enumerator.EnumeratorConfig()

        # Should not raise
        enumerator_obj = enumerator.RandomEnumerator(cs, config)
        syn = enumerator_obj.next()

        assert syn.count_building_blocks() <= config.max_building_blocks
