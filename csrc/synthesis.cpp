#include <algorithm>
#include <memory>
#include <sstream>

#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "molops.hpp"
#include "synthesis.hpp"
#include "utils/algorithm.hpp"
#include "utils/assert.hpp"
#include "utils/logging.hpp"

namespace synthesis_backend {

PostfixNotation::ItemType PostfixNotation::type(
    std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type index) const {
    Ensures(index < items.size());
    if (std::holds_alternative<Mol_sptr>(items[index])) {
        return ItemType::Molecule;
    } else {
        return ItemType::Reaction;
    }
}

void PostfixNotation::append(const Mol_sptr &mol) { items.push_back(mol); }

void PostfixNotation::append(const Reaction_sptr &reaction) {
    items.push_back(reaction);
}

void PostfixNotation::extend(const PostfixNotation &other) {
    items.insert(items.end(), other.items.begin(), other.items.end());
}

std::variant<Mol_sptr, Reaction_sptr> PostfixNotation::operator[](
    std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type index) const {
    Ensures(index < items.size());
    return items[index];
}

std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type
PostfixNotation::size() const {
    return items.size();
}

size_t PostfixNotation::count_reactions() const {
    return std::count_if(items.begin(), items.end(), [](const auto &item) {
        return std::holds_alternative<Reaction_sptr>(item);
    });
}

size_t PostfixNotation::count_building_blocks() const {
    return std::count_if(items.begin(), items.end(), [](const auto &item) {
        return std::holds_alternative<Mol_sptr>(item);
    });
}

std::string PostfixNotation::pickle() const {
    std::ostringstream buffer;
    boost::archive::binary_oarchive oa(buffer);
    size_t num_items = items.size();
    oa << num_items;
    for (const auto &item : items) {
        if (std::holds_alternative<Mol_sptr>(item)) {
            std::string item_pickle;
            RDKit::MolPickler::pickleMol(*std::get<Mol_sptr>(item), item_pickle,
                                         RDKit::PicklerOps::AllProps);
            oa << ItemType::Molecule << item_pickle;
        } else {
            std::string item_pickle;
            RDKit::ReactionPickler::pickleReaction(
                *std::get<Reaction_sptr>(item), item_pickle,
                RDKit::PicklerOps::AllProps);
            oa << ItemType::Reaction << item_pickle;
        }
    }
    return buffer.str();
}

std::unique_ptr<PostfixNotation> PostfixNotation::unpickle(const std::string &data) {
    std::istringstream buffer(data);
    boost::archive::binary_iarchive ia(buffer);
    size_t num_items;
    ia >> num_items;
    auto object = std::make_unique<PostfixNotation>();
    for (size_t i = 0; i < num_items; ++i) {
        PostfixNotation::ItemType type;
        std::string item_pickle;
        ia >> type >> item_pickle;
        if (type == ItemType::Molecule) {
            Mol_sptr mol(new RDKit::ROMol());
            RDKit::MolPickler::molFromPickle(item_pickle, mol.get(),
                                             RDKit::PicklerOps::AllProps);
            object->append(mol);
        } else {
            Reaction_sptr reaction =
                std::make_shared<Reaction_sptr::element_type>();
            RDKit::ReactionPickler::reactionFromPickle(item_pickle,
                                                       reaction.get());
            object->append(reaction);
        }
    }
    return object;
}

const PostfixNotation &Synthesis::get_postfix_notation() const {
    return postfix_notation;
}

const std::vector<MolSet> &Synthesis::get_stack() const { return stack; }

void Synthesis::push(const Mol_sptr &mol) {
    postfix_notation.append(mol);
    // Copy the molecule, because some RDKit operations might modify it,
    // although they are declared to be const.
    stack.push_back({Mol_sptr(new RDKit::ROMol(*mol))});
}

void Synthesis::push(const Reaction_sptr &reaction, size_t max_products) {
    if (stack.size() < reaction->getNumReactantTemplates()) {
        throw push_reaction_exception(
            "Not enough reactants for reaction. Expected " +
            std::to_string(reaction->getNumReactantTemplates()) +
            " reactants, but got " + std::to_string(stack.size()));
    }

    std::vector<MolSet> reactant_sets(
        stack.end() - reaction->getNumReactantTemplates(), stack.end());
    MolSet main_products;

    try {
        generate_combinations(
            reactant_sets, [&](std::vector<Mol_sptr> reactants) {
                std::sort(reactants.begin(), reactants.end());
                do {
                    auto product_groups = reaction->runReactants(reactants);
                    for (const auto &pgroup : product_groups) {
                        if (pgroup.empty()) {
                            continue;
                        }
                        auto product = sanitize(pgroup[0]);
                        if (product.has_value()) {
                            main_products.insert(product.value());
                            if (main_products.size() >= max_products) {
                                throw stop_iteration();
                            }
                        }
                    }
                } while (
                    std::next_permutation(reactants.begin(), reactants.end()));
            });
    } catch (const stop_iteration &) {
    }

    if (main_products.empty()) {
        throw push_reaction_exception(
            "No sanitized products generated from the reaction.");
    }

    postfix_notation.append(reaction);
    stack.erase(stack.end() - reaction->getNumReactantTemplates(), stack.end());
    stack.push_back(main_products);
}

void Synthesis::push(const Synthesis &synthesis) {
    postfix_notation.extend(synthesis.postfix_notation);
    stack.insert(stack.end(), synthesis.stack.begin(), synthesis.stack.end());
}

void Synthesis::push(const std::shared_ptr<Synthesis> &synthesis) {
    Ensures(synthesis != nullptr);
    push(*synthesis);
}

void Synthesis::push(
    const std::variant<Mol_sptr, std::shared_ptr<Synthesis>> &item) {
    if (std::holds_alternative<Mol_sptr>(item)) {
        push(std::get<Mol_sptr>(item));
    } else {
        push(std::get<std::shared_ptr<Synthesis>>(item));
    }
}

MolSet Synthesis::top(std::vector<MolSet>::size_type index) const {
    Ensures(index < stack.size());
    return stack[stack.size() - 1 - index];
}

MolSet Synthesis::top() const { return top(0); }

std::vector<MolSet>::size_type Synthesis::stack_size() const {
    return stack.size();
}

bool Synthesis::is_empty() const {
    return postfix_notation.size() == 0 && stack.empty();
}

size_t Synthesis::count_reactions() const {
    return postfix_notation.count_reactions();
}

size_t Synthesis::count_building_blocks() const {
    return postfix_notation.count_building_blocks();
}

std::string Synthesis::pickle() const {
    std::ostringstream buffer;
    boost::archive::binary_oarchive oa(buffer);
    std::string pfn_pickle = postfix_notation.pickle();
    oa << pfn_pickle;
    oa << stack.size();
    for (const auto &mol_set : stack) {
        size_t num_mols = mol_set.size();
        oa << num_mols;
        for (const auto &mol : mol_set) {
            std::string mol_pickle;
            RDKit::MolPickler::pickleMol(*mol, mol_pickle,
                                         RDKit::PicklerOps::AllProps);
            oa << mol_pickle;
        }
    }
    return buffer.str();
}

std::unique_ptr<Synthesis> Synthesis::unpickle(const std::string &data) {
    std::istringstream buffer(data);
    boost::archive::binary_iarchive ia(buffer);
    std::string pfn_pickle;
    ia >> pfn_pickle;
    size_t stack_size;
    ia >> stack_size;
    std::vector<MolSet> stack;
    for (size_t i = 0; i < stack_size; ++i) {
        size_t num_mols;
        ia >> num_mols;
        MolSet mol_set;
        for (size_t j = 0; j < num_mols; ++j) {
            std::string mol_pickle;
            ia >> mol_pickle;
            Mol_sptr mol(new RDKit::ROMol());
            RDKit::MolPickler::molFromPickle(mol_pickle, mol.get(),
                                             RDKit::PicklerOps::AllProps);
            mol_set.insert(mol);
        }
        stack.push_back(mol_set);
    }
    return std::make_unique<Synthesis>(std::move(*PostfixNotation::unpickle(pfn_pickle)),
                                        std::move(stack));
}

void save_synthesis_vector(const SynthesisVector &syntheses,
                           const std::filesystem::path &path) {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open file for writing");
    }
    boost::archive::binary_oarchive oa(ofs);
    size_t num_syntheses = syntheses.size();
    oa << num_syntheses;
    for (const auto &synthesis : syntheses) {
        std::string pickle = synthesis->pickle();
        oa << pickle;
    }
    ofs.close();
}

SynthesisVector *load_synthesis_vector(const std::filesystem::path &path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading");
    }
    boost::archive::binary_iarchive ia(ifs);
    size_t num_syntheses;
    ia >> num_syntheses;
    SynthesisVector *syntheses = new SynthesisVector();
    syntheses->reserve(num_syntheses);
    for (size_t i = 0; i < num_syntheses; ++i) {
        std::string pickle;
        ia >> pickle;
        std::unique_ptr<Synthesis> synthesis = Synthesis::unpickle(pickle);
        syntheses->emplace_back(std::move(synthesis));
    }
    ifs.close();
    return syntheses;
}
} // namespace synthesis_backend
