#pragma once

#include <functional>
#include <optional>
#include <vector>

#include "types.hpp"

namespace prexsyn_engine {

template <typename T>
std::vector<T> ecfp4_fingerprint(const std::optional<Mol_sptr> &mol);

template <typename T>
std::vector<T> rdkit_fingerprint(const std::optional<Mol_sptr> &mol);

template <typename T>
std::vector<T> fcfp4_fingerprint(const std::optional<Mol_sptr> &mol);

template <typename T>
using FpFunc = std::function<std::vector<T>(const std::optional<Mol_sptr> &)>;

template <typename T> FpFunc<T> get_fp_func(const std::string &name);

float tanimoto_similarity(const Mol_sptr &mol1, const Mol_sptr &mol2,
                          const std::string &fp_type);

} // namespace prexsyn_engine
