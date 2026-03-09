#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "postfix_notation.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

class ChemicalSpace;

class ChemicalSpaceSynthesis {
private:
    const ChemicalSpace &cs_;
    PostfixNotation postfix_notation_;
    std::shared_ptr<Synthesis> synthesis_;

public:
    ChemicalSpaceSynthesis(const ChemicalSpace &cs)
        : cs_(cs), postfix_notation_(), synthesis_(std::make_shared<Synthesis>()) {}

    const PostfixNotation &postfix_notation() const { return postfix_notation_; }
    const Synthesis &synthesis() const { return *synthesis_; }
    std::shared_ptr<Synthesis> &get_synthesis() { return synthesis_; }

    struct Result {
        bool is_ok;
        std::string message;

        static Result ok() { return {.is_ok = true, .message = {}}; }
        static Result error(const std::string &msg) { return {.is_ok = false, .message = msg}; }

        operator bool() const { return is_ok; }
    };

    size_t count_building_blocks() const;
    size_t count_reactions() const;
    std::vector<std::shared_ptr<Molecule>> products() const;

    Result add_building_block(BuildingBlockLibrary::Index) noexcept;
    Result add_building_block(const std::string &) noexcept;
    Result add_reaction(ReactionLibrary::Index) noexcept;
    Result add_reaction(const std::string &) noexcept;
    Result add_postfix_notation(const PostfixNotation &) noexcept;
    Result undo() noexcept;
};

using Synthesis = ChemicalSpaceSynthesis;

} // namespace prexsyn::chemspace
