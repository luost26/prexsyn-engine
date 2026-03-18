from pathlib import Path

import numpy as np

from prexsyn_engine import chemspace, datapipe, descriptor


def resource_path(name: str) -> Path:
    root = Path(__file__).resolve().parents[2]
    return root / "resources" / "test" / "chemspace_small_1" / name


def make_chemical_space():
    bb_lib = chemspace.bb_lib_from_sdf(resource_path("bb.sdf"))
    rxn_lib = chemspace.rxn_lib_from_plain_text(resource_path("rxn.txt"))
    int_lib = chemspace.IntermediateLibrary()

    cs = chemspace.ChemicalSpace(bb_lib, rxn_lib, int_lib)
    cs.build_reactant_lists_for_building_blocks()
    cs.generate_intermediates()
    cs.build_reactant_lists_for_intermediates()
    return cs


def test_generator_config_defaults_and_assignment():
    config = datapipe.GeneratorConfig()

    assert config.max_building_blocks == 5
    assert config.heavy_atom_limit == 50

    config.max_building_blocks = 3
    config.heavy_atom_limit = 24

    assert config.max_building_blocks == 3
    assert config.heavy_atom_limit == 24


def test_generator_next_returns_synthesis():
    cs = make_chemical_space()

    gen = datapipe.Generator(cs, random_seed=7)
    syn = gen.next()

    assert isinstance(syn, chemspace.Synthesis)
    assert syn.count_building_blocks() >= 1
    assert len(syn.products()) >= 1


def test_generator_next_with_product_returns_pair():
    cs = make_chemical_space()

    gen = datapipe.Generator(cs, random_seed=11)
    syn, product = gen.next_with_product()

    assert isinstance(syn, chemspace.Synthesis)
    assert product.smiles() != ""


def test_generator_same_seed_has_reproducible_first_product():
    cs = make_chemical_space()

    gen1 = datapipe.Generator(cs, random_seed=123)
    gen2 = datapipe.Generator(cs, random_seed=123)

    _, product1 = gen1.next_with_product()
    _, product2 = gen2.next_with_product()

    assert product1.smiles() == product2.smiles()


def test_data_pipeline_get_returns_expected_arrays():
    cs = make_chemical_space()

    molecule_descriptors = {"ecfp4": descriptor.MorganFingerprint.ecfp4()}
    synthesis_descriptors = {
        "pfn": descriptor.SynthesisPostfixNotation.create(max_length=8)
    }

    pipeline = datapipe.DataPipeline(cs, molecule_descriptors, synthesis_descriptors)
    pipeline.start_workers([42])
    try:
        batch = pipeline.get(2)
    finally:
        pipeline.stop_workers()

    # stop_workers should be safe even when workers are already stopped.
    pipeline.stop_workers()

    assert set(batch.keys()) == {"ecfp4", "pfn"}

    ecfp4 = batch["ecfp4"]
    pfn = batch["pfn"]

    assert isinstance(ecfp4, np.ndarray)
    assert ecfp4.dtype == np.dtype("bool")
    assert ecfp4.shape[0] == 2

    assert isinstance(pfn, np.ndarray)
    assert pfn.dtype == np.dtype("int64")
    assert pfn.shape == (2, 8, 3)


def test_two_data_pipeline_instances_can_run_concurrently():
    cs = make_chemical_space()

    pipeline1 = datapipe.DataPipeline(
        cs,
        {"ecfp4": descriptor.MorganFingerprint.ecfp4()},
        {"pfn": descriptor.SynthesisPostfixNotation.create(max_length=8)},
    )
    pipeline2 = datapipe.DataPipeline(
        cs,
        {"ecfp4": descriptor.MorganFingerprint.ecfp4()},
        {"pfn": descriptor.SynthesisPostfixNotation.create(max_length=8)},
    )

    pipeline1.start_workers([42])
    pipeline2.start_workers([43])
    try:
        batch1 = pipeline1.get(2)
        batch2 = pipeline2.get(2)
        print(batch1)
        print(batch2)
    finally:
        pipeline1.stop_workers()
        pipeline2.stop_workers()


def test_data_pipeline_stop_workers_without_start_is_safe():
    cs = make_chemical_space()

    pipeline = datapipe.DataPipeline(
        cs,
        {"ecfp4": descriptor.MorganFingerprint.ecfp4()},
        {"pfn": descriptor.SynthesisPostfixNotation.create(max_length=8)},
    )

    pipeline.stop_workers()
