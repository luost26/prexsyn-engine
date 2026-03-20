import argparse
from pathlib import Path

import prexsyn_engine

_default_sdf_path = (
    "data/Enamine_Rush-Delivery_Building_Blocks-US_223244cmpd_20231001.sdf"
)
_default_bb_lib_cache = "data/enamine.bin"
_default_rxn_path = "data/reactions/rxn115.txt"
_default_rxn_cache = "data/rxn115.bin"
_default_output_path = "data/enamine_rxn115_chemspace.bin"


def load_building_blocks(sdf_path=_default_sdf_path, cache_path=_default_bb_lib_cache):
    cache_path = Path(cache_path)
    if cache_path.exists():
        print(f"Loading building block library from cache at {cache_path}")
        bb_lib = prexsyn_engine.chemspace.BuildingBlockLibrary.deserialize(cache_path)
    else:
        bb_lib = prexsyn_engine.chemspace.bb_lib_from_sdf(sdf_path)
        bb_lib.serialize(cache_path)

    print(f"Loaded {len(bb_lib)} building blocks")
    return bb_lib


def load_reactions(rxn_path=_default_rxn_path, cache_path=_default_rxn_cache):
    cache_path = Path(cache_path)
    if cache_path.exists():
        print(f"Loading reaction library from cache at {cache_path}")
        rxn_lib = prexsyn_engine.chemspace.ReactionLibrary.deserialize(cache_path)
    else:
        rxn_lib = prexsyn_engine.chemspace.rxn_lib_from_plain_text(rxn_path)
        rxn_lib.serialize(cache_path)

    print(f"Loaded {len(rxn_lib)} reactions")
    return rxn_lib


def create_empty_intermediate_library():
    int_lib = prexsyn_engine.chemspace.IntermediateLibrary()
    print(f"Created empty intermediate library with {len(int_lib)} intermediates")
    return int_lib


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build chemical space for Enamine-Rxn115 dataset"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rebuild of chemical space even if output file exists",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=_default_output_path,
        help=f"Path to save the chemical space (default: {_default_output_path})",
    )
    parser.add_argument(
        "--sdf",
        type=str,
        default=_default_sdf_path,
        help=f"Path to building block SDF file (default: {_default_sdf_path})",
    )
    parser.add_argument(
        "--rxn",
        type=str,
        default=_default_rxn_path,
        help=f"Path to reaction text file (default: {_default_rxn_path})",
    )

    args = parser.parse_args()

    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        peek_stats = prexsyn_engine.chemspace.ChemicalSpace.peek(output_path)
        print(f"Chemical space at {output_path} already exists:")
        print("- Number of building blocks:", peek_stats.num_building_blocks)
        print("- Number of reactions:", peek_stats.num_reactions)
        print("- Number of intermediates:", peek_stats.num_intermediates)
    else:
        bb_lib = load_building_blocks(sdf_path=args.sdf)
        rxn_lib = load_reactions(rxn_path=args.rxn)
        int_lib = create_empty_intermediate_library()

        cs = prexsyn_engine.chemspace.ChemicalSpace(
            bb_lib=bb_lib, rxn_lib=rxn_lib, int_lib=int_lib
        )
        cs.build_reactant_lists_for_building_blocks()
        cs.generate_intermediates()
        cs.build_reactant_lists_for_intermediates()
        cs.serialize(output_path)
