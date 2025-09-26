#pragma once

#include <memory>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "../types.hpp"

namespace synthesis_backend {
inline std::shared_ptr<spdlog::logger> logger() {
    static auto logger = spdlog::get("synthesis_backend");
    if (!logger) {
        logger = spdlog::stderr_color_mt("synthesis_backend");
    }
    return logger;
}
} // namespace synthesis_backend

template <>
struct fmt::formatter<synthesis_backend::Mol_sptr>
    : fmt::formatter<std::string> {
    auto format(synthesis_backend::Mol_sptr mol, format_context &ctx) const
        -> decltype(ctx.out()) {
        if (!mol) {
            return fmt::format_to(ctx.out(), "Mol(null)");
        } else {
            return fmt::format_to(ctx.out(), "Mol({})",
                                  RDKit::MolToSmiles(*mol));
        }
    }
};

template <>
struct fmt::formatter<synthesis_backend::Reaction_sptr>
    : fmt::formatter<std::string> {
    auto format(synthesis_backend::Reaction_sptr reaction,
                format_context &ctx) const -> decltype(ctx.out()) {
        if (!reaction) {
            return fmt::format_to(ctx.out(), "Reaction(null)");
        } else {
            return fmt::format_to(
                ctx.out(), "Reaction({})",
                RDKit::ChemicalReactionToRxnSmarts(*reaction));
        }
    }
};
