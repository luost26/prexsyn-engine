#include "synthesis.hpp"

#include <cstddef>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "chemical_space.hpp"
#include "int_lib.hpp"
#include "postfix_notation.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

size_t ChemicalSpaceSynthesis::count_building_blocks() const {
    size_t count = 0;
    for (const auto &token : postfix_notation_.tokens()) {
        if (token.type == PostfixNotation::Token::Type::BuildingBlock) {
            ++count;
        }
    }
    return count;
}

size_t ChemicalSpaceSynthesis::count_reactions() const {
    size_t count = 0;
    for (const auto &token : postfix_notation_.tokens()) {
        if (token.type == PostfixNotation::Token::Type::Reaction) {
            ++count;
        }
    }
    return count;
}

std::vector<std::shared_ptr<Molecule>> ChemicalSpaceSynthesis::products() const {
    const auto &top = synthesis_->stack_top();
    std::vector<std::shared_ptr<Molecule>> result;
    result.reserve(top->size());
    for (size_t i = 0; i < top->size(); ++i) {
        result.push_back(top->at(i));
    }
    return result;
}

using Result = ChemicalSpaceSynthesis::Result;

Result ChemicalSpaceSynthesis::add_building_block(BuildingBlockLibrary::Index index) noexcept {
    try {
        const auto &bb_item = cs_.bb_lib().get(index);
        synthesis_->push(bb_item.molecule);
        postfix_notation_.append(bb_item.index, PostfixNotation::Token::Type::BuildingBlock);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

Result ChemicalSpaceSynthesis::add_building_block(const std::string &index) noexcept {
    try {
        const auto &bb_item = cs_.bb_lib().get(index);
        synthesis_->push(bb_item.molecule);
        postfix_notation_.append(bb_item.index, PostfixNotation::Token::Type::BuildingBlock);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

Result ChemicalSpaceSynthesis::add_reaction(ReactionLibrary::Index index) noexcept {
    try {
        const auto &rxn_item = cs_.rxn_lib().get(index);
        synthesis_->push(rxn_item.reaction);
        postfix_notation_.append(rxn_item.index, PostfixNotation::Token::Type::Reaction);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

Result ChemicalSpaceSynthesis::add_reaction(const std::string &index) noexcept {
    try {
        const auto &rxn_item = cs_.rxn_lib().get(index);
        synthesis_->push(rxn_item.reaction);
        postfix_notation_.append(rxn_item.index, PostfixNotation::Token::Type::Reaction);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

Result ChemicalSpaceSynthesis::add_postfix_notation(const PostfixNotation &pfn) noexcept {
    size_t count_operations = 0;
    auto result = Result::ok();
    for (const auto &token : pfn.tokens()) {
        if (token.type == PostfixNotation::Token::Type::BuildingBlock) {
            result = add_building_block(token.index);
            if (result) {
                count_operations++;
            } else {
                break;
            }
        } else if (token.type == PostfixNotation::Token::Type::Reaction) {
            result = add_reaction(token.index);
            if (result) {
                count_operations++;
            } else {
                break;
            }
        }
    }

    if (!result) {
        for (size_t i = 0; i < count_operations; ++i) {
            undo();
        }
    }
    return result;
};

Result ChemicalSpaceSynthesis::add_intermediate(IntermediateLibrary::Index index) noexcept {
    try {
        const auto &int_item = cs_.int_lib().get(index);
        return add_postfix_notation(int_item.postfix_notation);
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

Result ChemicalSpaceSynthesis::undo() noexcept {
    try {
        postfix_notation_.pop_back();
        synthesis_->undo();
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

} // namespace prexsyn::chemspace
