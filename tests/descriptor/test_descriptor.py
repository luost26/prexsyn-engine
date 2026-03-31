import numpy as np
import pytest

from prexsyn_engine import chemistry, descriptor


Molecule = chemistry.Molecule


class TestMorganFingerprint:
    """Test cases for MorganFingerprint descriptor."""

    def test_ecfp4_creation_succeeds(self):
        """Test that ECFP4 fingerprint can be created."""
        fp = descriptor.MorganFingerprint.ecfp4()
        assert fp is not None

    def test_ecfp4_dtype_is_bool8(self):
        """Test that ECFP4 returns bool8 dtype."""
        fp = descriptor.MorganFingerprint.ecfp4()
        assert fp.dtype() == np.dtype("bool")

    def test_ecfp4_size_returns_tuple(self):
        """Test that size() returns a tuple."""
        fp = descriptor.MorganFingerprint.ecfp4()
        size = fp.size()
        assert isinstance(size, (tuple, list))
        assert len(size) > 0
        assert isinstance(size[0], int)

    def test_ecfp4_num_elements_matches_size(self):
        """Test that num_elements() matches product of size()."""
        fp = descriptor.MorganFingerprint.ecfp4()
        size = fp.size()
        expected_num_elements = 1
        for s in size:
            expected_num_elements *= s
        assert fp.num_elements() == expected_num_elements

    def test_ecfp4_size_in_bytes_is_positive(self):
        """Test that size_in_bytes() returns positive value."""
        fp = descriptor.MorganFingerprint.ecfp4()
        assert fp.size_in_bytes() > 0

    def test_ecfp4_call_returns_numpy_array(self):
        """Test that calling ECFP4 on a molecule returns numpy array."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCO")

        result = fp(mol)

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.dtype("bool")

    def test_ecfp4_call_array_size_matches_descriptor(self):
        """Test that returned array size matches descriptor size."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCO")

        result = fp(mol)
        expected_shape = fp.size()

        assert result.shape == tuple(expected_shape)

    def test_ecfp4_different_molecules_produce_different_fingerprints(self):
        """Test that different molecules produce different ECFP4 fingerprints."""
        fp = descriptor.MorganFingerprint.ecfp4()

        mol1_fp = fp(Molecule.from_smiles("CCO"))
        mol2_fp = fp(Molecule.from_smiles("CC"))

        # At least one bit should differ
        assert not np.array_equal(mol1_fp, mol2_fp)

    def test_ecfp4_same_molecule_produces_same_fingerprint(self):
        """Test that the same molecule produces the same ECFP4 fingerprint."""
        fp = descriptor.MorganFingerprint.ecfp4()
        smiles = "CCO"

        fp1 = fp(Molecule.from_smiles(smiles))
        fp2 = fp(Molecule.from_smiles(smiles))

        assert np.array_equal(fp1, fp2)

    def test_fcfp4_creation_succeeds(self):
        """Test that FCFP4 fingerprint can be created."""
        fp = descriptor.MorganFingerprint.fcfp4()
        assert fp is not None


class TestTokenDef:
    """Test cases for TokenDef class."""

    def test_token_def_creation_succeeds(self):
        """Test that TokenDef can be created."""
        token_def = descriptor.TokenDef()
        assert token_def is not None

    def test_token_def_default_values(self):
        """Test that TokenDef has expected default values."""
        token_def = descriptor.TokenDef()

        assert token_def.pad == 0
        assert token_def.end == 1
        assert token_def.start == 2
        assert token_def.bb == 3
        assert token_def.rxn == 4

    def test_token_def_attributes_are_writeable(self):
        """Test that TokenDef attributes can be modified."""
        token_def = descriptor.TokenDef()

        token_def.pad = 10
        token_def.end = 11
        token_def.start = 12
        token_def.bb = 13
        token_def.rxn = 14

        assert token_def.pad == 10
        assert token_def.end == 11
        assert token_def.start == 12
        assert token_def.bb == 13
        assert token_def.rxn == 14

    def test_token_def_from_dict_with_all_keys(self):
        """Test creating TokenDef from dict with all keys."""
        d = {"pad": 5, "end": 6, "start": 7, "bb": 8, "rxn": 9}

        token_def = descriptor.TokenDef.from_dict(d)

        assert token_def.pad == 5
        assert token_def.end == 6
        assert token_def.start == 7
        assert token_def.bb == 8
        assert token_def.rxn == 9

    def test_token_def_from_dict_with_partial_keys(self):
        """Test creating TokenDef from dict with partial keys."""
        d = {"pad": 5, "bb": 8}

        token_def = descriptor.TokenDef.from_dict(d)

        assert token_def.pad == 5
        assert token_def.end == 1  # default
        assert token_def.start == 2  # default
        assert token_def.bb == 8
        assert token_def.rxn == 4  # default

    def test_token_def_from_dict_with_empty_dict(self):
        """Test creating TokenDef from empty dict."""
        token_def = descriptor.TokenDef.from_dict({})

        assert token_def.pad == 0
        assert token_def.end == 1
        assert token_def.start == 2
        assert token_def.bb == 3
        assert token_def.rxn == 4

    def test_token_def_from_dict_with_zero_values(self):
        """Test creating TokenDef from dict with zero values."""
        d = {"pad": 0, "end": 0, "start": 0, "bb": 0, "rxn": 0}

        token_def = descriptor.TokenDef.from_dict(d)

        assert token_def.pad == 0
        assert token_def.end == 0
        assert token_def.start == 0
        assert token_def.bb == 0
        assert token_def.rxn == 0


class TestSynthesisPostfixNotation:
    """Test cases for SynthesisPostfixNotation descriptor."""

    def test_create_with_default_args_succeeds(self):
        """Test creating SynthesisPostfixNotation with default arguments."""
        desc = descriptor.SynthesisPostfixNotation.create()
        assert desc is not None

    def test_create_with_max_length_succeeds(self):
        """Test creating SynthesisPostfixNotation with max_length."""
        desc = descriptor.SynthesisPostfixNotation.create(max_length=32)
        assert desc is not None

    def test_create_with_token_def_succeeds(self):
        """Test creating SynthesisPostfixNotation with TokenDef."""
        token_def = descriptor.TokenDef()
        desc = descriptor.SynthesisPostfixNotation.create(token_def)
        assert desc is not None

    def test_create_with_token_def_and_max_length_succeeds(self):
        """Test creating SynthesisPostfixNotation with TokenDef and max_length."""
        token_def = descriptor.TokenDef()
        desc = descriptor.SynthesisPostfixNotation.create(token_def, max_length=32)
        assert desc is not None

    def test_dtype_is_int64(self):
        """Test that dtype() returns int64."""
        desc = descriptor.SynthesisPostfixNotation.create()
        assert desc.dtype() == np.dtype("int64")

    def test_size_returns_tuple(self):
        """Test that size() returns a tuple."""
        desc = descriptor.SynthesisPostfixNotation.create()
        size = desc.size()
        assert isinstance(size, (tuple, list))
        assert len(size) == 2

    def test_size_shape_for_default_max_length(self):
        """Test size() shape for default max_length of 16."""
        desc = descriptor.SynthesisPostfixNotation.create()
        size = desc.size()

        assert size[0] == 16  # default max_length
        assert size[1] == 3  # type, bb_idx, rxn_idx

    def test_size_shape_for_custom_max_length(self):
        """Test size() shape for custom max_length."""
        max_length = 32
        desc = descriptor.SynthesisPostfixNotation.create(max_length=max_length)
        size = desc.size()

        assert size[0] == max_length
        assert size[1] == 3

    def test_num_elements_matches_size(self):
        """Test that num_elements() matches product of size()."""
        desc = descriptor.SynthesisPostfixNotation.create()
        size = desc.size()
        expected = size[0] * size[1]

        assert desc.num_elements() == expected

    def test_size_in_bytes_is_positive(self):
        """Test that size_in_bytes() is positive."""
        desc = descriptor.SynthesisPostfixNotation.create()
        # Each int64 is 8 bytes
        expected_bytes = desc.num_elements() * 8
        assert desc.size_in_bytes() == expected_bytes

    def test_size_in_bytes_increases_with_max_length(self):
        """Test that size_in_bytes() increases with max_length."""
        desc_small = descriptor.SynthesisPostfixNotation.create(max_length=8)
        desc_large = descriptor.SynthesisPostfixNotation.create(max_length=32)

        assert desc_large.size_in_bytes() > desc_small.size_in_bytes()

    def test_min_max_length_is_4(self):
        """Test that max_length < 4 raises ValueError."""
        with pytest.raises((ValueError, RuntimeError)):
            descriptor.SynthesisPostfixNotation.create(max_length=3)

    def test_max_length_4_succeeds(self):
        """Test that max_length of 4 is accepted."""
        desc = descriptor.SynthesisPostfixNotation.create(max_length=4)
        assert desc.size()[0] == 4

    def test_custom_token_def_is_used(self):
        """Test that custom TokenDef is used in descriptor."""
        # Just verify that we can create with custom token def
        token_def = descriptor.TokenDef()
        token_def.pad = 100
        token_def.end = 101

        desc = descriptor.SynthesisPostfixNotation.create(token_def)
        assert desc is not None

    def test_create_overload_with_token_def_only(self):
        """Test create() with only TokenDef argument."""
        token_def = descriptor.TokenDef()
        desc = descriptor.SynthesisPostfixNotation.create(token_def)

        # Should use default max_length of 16
        assert desc.size()[0] == 16

    def test_create_multiple_instances_are_independent(self):
        """Test that multiple instances are independent."""
        desc1 = descriptor.SynthesisPostfixNotation.create(max_length=16)
        desc2 = descriptor.SynthesisPostfixNotation.create(max_length=32)

        assert desc1.size()[0] == 16
        assert desc2.size()[0] == 32


class TestBatchedMorganFingerprintCalls:
    """Test cases for batched __call__ on MorganFingerprint descriptor."""

    def test_batched_ecfp4_with_single_molecule_returns_2d_array(self):
        """Test that calling ECFP4 with single molecule in list returns 2D array."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCC")

        result = fp([mol])

        assert isinstance(result, np.ndarray)
        assert len(result.shape) == 2
        assert result.shape[0] == 1

    def test_batched_ecfp4_with_multiple_molecules_returns_correct_shape(self):
        """Test that batched call returns shape (batch_size, *descriptor_size)."""
        fp = descriptor.MorganFingerprint.ecfp4()
        molecules = [
            Molecule.from_smiles("CCC"),
            Molecule.from_smiles("CCO"),
            Molecule.from_smiles("CCCC"),
        ]

        result = fp(molecules)

        expected_shape = (len(molecules),) + tuple(fp.size())
        assert result.shape == expected_shape

    def test_batched_ecfp4_dtype_matches_single_call(self):
        """Test that batched call dtype matches single call."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCC")

        single_result = fp(mol)
        batched_result = fp([mol])

        assert batched_result.dtype == single_result.dtype

    def test_batched_ecfp4_individual_results_match_single_calls(self):
        """Test that each batch result matches individual single calls."""
        fp = descriptor.MorganFingerprint.ecfp4()
        smiles_list = ["CCC", "CCO", "CCCC", "CC"]
        molecules = [Molecule.from_smiles(s) for s in smiles_list]

        # Get individual results
        individual_results = [fp(mol) for mol in molecules]

        # Get batched result
        batched_result = fp(molecules)

        # Each batch element should match the corresponding individual result
        for i, individual in enumerate(individual_results):
            assert np.array_equal(batched_result[i], individual), (
                f"Batch result at index {i} does not match individual result"
            )

    def test_batched_fcfp4_individual_results_match_single_calls(self):
        """Test batched FCFP4 results match individual calls."""
        fp = descriptor.MorganFingerprint.fcfp4()
        molecules = [
            Molecule.from_smiles("CCC"),
            Molecule.from_smiles("CCO"),
            Molecule.from_smiles("CC"),
        ]

        individual_results = [fp(mol) for mol in molecules]
        batched_result = fp(molecules)

        for i, individual in enumerate(individual_results):
            assert np.array_equal(batched_result[i], individual)

    def test_batched_ecfp4_empty_list_returns_correct_shape(self):
        """Test that empty molecule list returns (0, *fingerprint_size)."""
        fp = descriptor.MorganFingerprint.ecfp4()

        result = fp([])

        assert result.shape[0] == 0
        assert len(result.shape) == len(fp.size()) + 1

    def test_batched_ecfp4_large_batch_consistency(self):
        """Test batched call consistency with larger batch."""
        fp = descriptor.MorganFingerprint.ecfp4()
        molecules = [Molecule.from_smiles("CCC") for _ in range(10)]

        individual_results = [fp(mol) for mol in molecules]
        batched_result = fp(molecules)

        for i, individual in enumerate(individual_results):
            assert np.array_equal(batched_result[i], individual)

    def test_batched_call_all_molecules_produce_arrays(self):
        """Test that all batch results are valid numpy arrays."""
        fp = descriptor.MorganFingerprint.ecfp4()
        molecules = [
            Molecule.from_smiles("CCC"),
            Molecule.from_smiles("c1ccccc1"),  # benzene
            Molecule.from_smiles("CCN"),
        ]

        result = fp(molecules)

        assert result.dtype == np.dtype("bool")
        assert result.shape[0] == 3
        # Each result should be a valid array without NaN or inf
        assert np.all(np.isfinite(result.astype(int)))

    def test_batched_ecfp4_vs_single_produces_same_values(self):
        """Test that batched ECFP4 produces identical values vs calling individually."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCO")

        # Call as batch of 1
        batched_single = fp([mol])

        # Call individually
        individual = fp(mol)

        # Reshape individual to match batch
        individual_as_batch = individual[np.newaxis, :]

        assert np.array_equal(batched_single, individual_as_batch)

    def test_batched_different_molecules_produce_different_fingerprints(self):
        """Test that different molecules in batch produce different fingerprints."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol1 = Molecule.from_smiles("CCO")
        mol2 = Molecule.from_smiles("CC")

        result = fp([mol1, mol2])

        # Results should have different fingerprints
        assert not np.array_equal(result[0], result[1])

    def test_batched_identical_molecules_produce_identical_fingerprints(self):
        """Test that identical molecules in batch produce identical fingerprints."""
        fp = descriptor.MorganFingerprint.ecfp4()
        mol = Molecule.from_smiles("CCCC")

        result = fp([mol, mol, mol])

        # All should be identical
        assert np.array_equal(result[0], result[1])
        assert np.array_equal(result[1], result[2])
