#pragma once

#include <memory>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "../types.hpp"

namespace prexsyn_engine {
inline std::shared_ptr<spdlog::logger> logger() {
    auto logger = spdlog::get("prexsyn_engine");
    if (!logger) {
        logger = spdlog::stderr_color_mt("prexsyn_engine");
    }
    return logger;
}
} // namespace prexsyn_engine

template <>
struct fmt::formatter<prexsyn_engine::Mol_sptr> : fmt::formatter<std::string> {
    auto format(prexsyn_engine::Mol_sptr mol, format_context &ctx) const
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
struct fmt::formatter<prexsyn_engine::Reaction_sptr>
    : fmt::formatter<std::string> {
    auto format(prexsyn_engine::Reaction_sptr reaction,
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
