from pathlib import Path

import numpy as np
import pytest

from prexsyn_engine import chemspace, descriptor, detokenizer


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


def make_valid_tokens(
    cs: chemspace.ChemicalSpace, token_def: descriptor.TokenDef
) -> np.ndarray:
    bb_idx_a = cs.bb_lib().get("EN300-250786").index
    bb_idx_b = cs.bb_lib().get("EN300-101318").index
    rxn_idx = cs.rxn_lib().get("ReactionA").index

    return np.array(
        [
            [token_def.bb, bb_idx_a, 0],
            [token_def.bb, bb_idx_b, 0],
            [token_def.rxn, 0, rxn_idx],
            [token_def.end, 0, 0],
        ],
        dtype=np.int64,
    )


def extract_events(tokens: np.ndarray, token_def: descriptor.TokenDef):
    events = []
    for row in tokens:
        token_type = int(row[0])
        if token_type in (token_def.pad, token_def.start):
            continue
        if token_type == token_def.end:
            break
        if token_type == token_def.bb:
            events.append(("bb", int(row[1])))
        elif token_type == token_def.rxn:
            events.append(("rxn", int(row[2])))
    return events


def make_valid_synthesis(cs: chemspace.ChemicalSpace):
    syn = cs.new_synthesis()
    assert syn.add_building_block("EN300-250786").is_ok
    assert syn.add_building_block("EN300-101318").is_ok
    assert syn.add_reaction("ReactionA", None).is_ok
    return syn


def test_detokenize_single_sequence_returns_valid_synthesis():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()
    tokens = make_valid_tokens(cs, token_def)

    syn = detokenizer.detokenize(tokens, cs, token_def)

    assert syn.count_building_blocks() == 2
    assert syn.count_reactions() == 1
    assert len(syn.products()) > 0


def test_detokenize_uses_default_token_def_when_omitted():
    cs = make_chemical_space()
    tokens = make_valid_tokens(cs, descriptor.TokenDef())

    syn = detokenizer.detokenize(tokens, cs)

    assert syn.count_building_blocks() == 2
    assert syn.count_reactions() == 1


def test_detokenize_raises_for_invalid_shape():
    cs = make_chemical_space()
    bad_tokens = np.array([1, 2, 3], dtype=np.int64)

    with pytest.raises((ValueError, RuntimeError), match="2D array"):
        detokenizer.detokenize(bad_tokens, cs)


def test_multithreaded_detokenizer_batch_decodes_each_item():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()

    tokens = make_valid_tokens(cs, token_def)
    batch_tokens = np.stack([tokens, tokens], axis=0)

    mt_detok = detokenizer.MultiThreadedDetokenizer(cs, token_def)
    outputs = mt_detok(2, batch_tokens)

    assert len(outputs) == 2
    for syn in outputs:
        assert syn.count_building_blocks() == 2
        assert syn.count_reactions() == 1
        assert len(syn.products()) > 0


def test_multithreaded_detokenizer_raises_for_invalid_shape():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()

    bad_batch = np.zeros((2, 4), dtype=np.int64)
    mt_detok = detokenizer.MultiThreadedDetokenizer(cs, token_def)

    with pytest.raises((ValueError, RuntimeError), match="3D array"):
        mt_detok(2, bad_batch)


def test_detokenize_then_tokenize_round_trip_is_consistent():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()
    input_tokens = make_valid_tokens(cs, token_def)

    syn = detokenizer.detokenize(input_tokens, cs, token_def)
    tokenizer = descriptor.SynthesisPostfixNotation.create(token_def, max_length=8)
    output_tokens = tokenizer(syn)

    assert output_tokens.shape == (8, 3)
    assert int(output_tokens[0, 0]) == token_def.start
    assert extract_events(output_tokens, token_def) == extract_events(
        input_tokens, token_def
    )


def test_multithreaded_detokenize_then_tokenize_round_trip_is_consistent():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()
    input_tokens = make_valid_tokens(cs, token_def)
    batch_tokens = np.stack([input_tokens, input_tokens], axis=0)

    mt_detok = detokenizer.MultiThreadedDetokenizer(cs, token_def)
    syntheses = mt_detok(2, batch_tokens)
    tokenizer = descriptor.SynthesisPostfixNotation.create(token_def, max_length=8)

    for syn in syntheses:
        output_tokens = tokenizer(syn)
        assert output_tokens.shape == (8, 3)
        assert int(output_tokens[0, 0]) == token_def.start
        assert extract_events(output_tokens, token_def) == extract_events(
            input_tokens, token_def
        )


def test_tokenize_then_detokenize_round_trip_is_consistent():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()
    tokenizer = descriptor.SynthesisPostfixNotation.create(token_def, max_length=8)

    input_synthesis = make_valid_synthesis(cs)
    input_tokens = tokenizer(input_synthesis)

    output_synthesis = detokenizer.detokenize(input_tokens, cs, token_def)
    output_tokens = tokenizer(output_synthesis)

    assert output_tokens.shape == (8, 3)
    assert extract_events(output_tokens, token_def) == extract_events(
        input_tokens, token_def
    )
    assert len(output_synthesis.products()) > 0


def test_tokenize_then_multithreaded_detokenize_round_trip_is_consistent():
    cs = make_chemical_space()
    token_def = descriptor.TokenDef()
    tokenizer = descriptor.SynthesisPostfixNotation.create(token_def, max_length=8)

    input_synthesis = make_valid_synthesis(cs)
    input_tokens = tokenizer(input_synthesis)
    batch_tokens = np.stack([input_tokens, input_tokens], axis=0)

    mt_detok = detokenizer.MultiThreadedDetokenizer(cs, token_def)
    output_syntheses = mt_detok(2, batch_tokens)

    assert len(output_syntheses) == 2
    for syn in output_syntheses:
        output_tokens = tokenizer(syn)
        assert output_tokens.shape == (8, 3)
        assert extract_events(output_tokens, token_def) == extract_events(
            input_tokens, token_def
        )
        assert len(syn.products()) > 0
