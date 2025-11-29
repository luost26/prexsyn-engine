import pathlib

from prexsyn_engine import chemspace, featurizer, pipeline

if __name__ == "__main__":
    cache_path = pathlib.Path("data/csd_cache_1")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if not cache_path.exists():
        csd = (
            chemspace.ChemicalSpaceDefinitionBuilder()
            .building_blocks_from_sdf("resources/test/chemspace_small_1/bb.sdf")
            .reactions_from_txt("resources/test/chemspace_small_1/rxn.txt")
            .secondary_building_blocks_from_single_reaction()
            .build_primary_index()
            .build_secondary_index()
            .build()
        )
        csd.save(str(cache_path))
    else:
        csd = chemspace.ChemicalSpaceDefinitionBuilder().all_from_cache(str(cache_path)).build()

    for rxn in csd.get_reactions():
        print(rxn.GetProp("original_index"), rxn.GetProp("reaction_index"))

    fs = (
        featurizer.FeaturizerSet()
        .add(featurizer.fingerprint.FingerprintFeaturizer("ecfp4", "ecfp4"))
        .add(featurizer.fingerprint.FingerprintFeaturizer("fcfp4", "fcfp4"))
        .add(featurizer.scaffold.MurckoScaffoldFeaturizer("murcko_scaffold", "fcfp4"))
        .add(featurizer.substructures.BRICSFragmentsFeaturizer("brics", "ecfp4", 8))
        .add(featurizer.synthesis.PostfixNotationFeaturizer())
    )

    gen_option = chemspace.SynthesisGeneratorOption()

    ppl = pipeline.DataPipelineV2(
        num_threads=16,
        csd=csd,
        gen_option=gen_option,
        featurizer=fs,
    )
    ppl.start()

    for i in range(1_000_000):
        data = ppl.get(4)
        print(data)
        print(data["ecfp4.fingerprint"].sum(-1))
        input("Press enter ...")
