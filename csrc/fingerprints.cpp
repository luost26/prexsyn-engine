#include "fingerprints.hpp"

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <memory>

#include "utils/assert.hpp"

namespace prexsyn_engine {

template <typename T>
std::vector<T> ecfp4_fingerprint(const std::optional<Mol_sptr> &mol) {
    std::vector<T> vec(2048);
    if (!mol.has_value()) {
        return vec;
    }
    Ensures(mol.value() != nullptr);
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> fp_gen(
        RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(
            2, false, false, true, false, true, nullptr, nullptr));
    // Explicitly set FFA is to avoid memory leakage in some RDKit versions
    // because if FFA is not specified, the function call would match to an
    // overload which returns a raw pointer to bit vector
    RDKit::FingerprintFuncArguments ffa{nullptr, nullptr, -1,
                                        nullptr, nullptr, nullptr};
    auto fp = fp_gen->getFingerprint(*mol.value(), ffa);
    Ensures(fp->size() == vec.size());
    for (size_t i = 0; i < fp->size(); ++i) {
        vec[i] = (*fp)[i];
    }
    return vec;
}

template std::vector<float>
ecfp4_fingerprint<float>(const std::optional<Mol_sptr> &mol);

template <typename T>
std::vector<T> rdkit_fingerprint(const std::optional<Mol_sptr> &mol) {
    std::vector<T> vec(2048);
    if (!mol.has_value()) {
        return vec;
    }
    Ensures(mol.value() != nullptr);
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> fp_gen(
        RDKit::RDKitFP::getRDKitFPGenerator<std::uint64_t>(
            1, 7, true, true, true, nullptr, false, {1, 2, 4, 8}, 2048, 2,
            false));
    RDKit::FingerprintFuncArguments ffa{nullptr, nullptr, -1,
                                        nullptr, nullptr, nullptr};
    auto fp = fp_gen->getFingerprint(*mol.value(), ffa);
    Ensures(fp->size() == vec.size());
    for (size_t i = 0; i < fp->size(); ++i) {
        vec[i] = (*fp)[i];
    }
    return vec;
}

template std::vector<float>
rdkit_fingerprint<float>(const std::optional<Mol_sptr> &mol);

template <typename T>
std::vector<T> fcfp4_fingerprint(const std::optional<Mol_sptr> &mol) {
    std::vector<T> vec(2048);
    if (!mol.has_value()) {
        return vec;
    }
    Ensures(mol.value() != nullptr);
    auto feature_inv_gen =
        RDKit::MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> fp_gen(
        RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(
            2, false, false, true, false, true, &feature_inv_gen, nullptr));
    // Explicitly set FFA is to avoid memory leakage in some RDKit versions
    // because if FFA is not specified, the function call would match to an
    // overload which returns a raw pointer to bit vector
    RDKit::FingerprintFuncArguments ffa{nullptr, nullptr, -1,
                                        nullptr, nullptr, nullptr};
    auto fp = fp_gen->getFingerprint(*mol.value(), ffa);
    Ensures(fp->size() == vec.size());
    for (size_t i = 0; i < fp->size(); ++i) {
        vec[i] = (*fp)[i];
    }
    return vec;
}

template std::vector<float>
fcfp4_fingerprint<float>(const std::optional<Mol_sptr> &mol);

template <typename T> FpFunc<T> get_fp_func(const std::string &name) {
    if (name == "ecfp4") {
        return ecfp4_fingerprint<T>;
    } else if (name == "rdkit") {
        return rdkit_fingerprint<T>;
    } else if (name == "fcfp4") {
        return fcfp4_fingerprint<T>;
    }
    throw std::runtime_error("Unknown fingerprint type: " + name);
}

template FpFunc<float> get_fp_func<float>(const std::string &name);

float tanimoto_similarity(const Mol_sptr &mol1, const Mol_sptr &mol2,
                          const std::string &fp_type) {
    Ensures(mol1 != nullptr && mol2 != nullptr);
    auto fp1 = get_fp_func<float>(fp_type)(mol1);
    auto fp2 = get_fp_func<float>(fp_type)(mol2);
    Ensures(fp1.size() == fp2.size());

    float intersection = 0;
    float union_size = 0;
    for (size_t i = 0; i < fp1.size(); ++i) {
        intersection += std::min(fp1[i], fp2[i]);
        union_size += std::max(fp1[i], fp2[i]);
    }
    if (union_size == 0.0f) {
        return 0.0f;
    }
    return intersection / union_size;
}

} // namespace prexsyn_engine
