#pragma once

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <variant>
#include <vector>

#include <GraphMol/GraphMol.h>

#include "types.hpp"

namespace prexsyn_engine {
class PostfixNotation {
    std::vector<std::variant<Mol_sptr, Reaction_sptr>> items;

  public:
    enum class ItemType { Molecule, Reaction };
    ItemType type(
        std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type) const;
    void append(const Mol_sptr &mol);
    void append(const Reaction_sptr &reaction);
    void extend(const PostfixNotation &other);
    std::variant<Mol_sptr, Reaction_sptr> operator[](
        std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type index)
        const;
    std::vector<std::variant<Mol_sptr, Reaction_sptr>>::size_type size() const;
    size_t count_reactions() const;
    size_t count_building_blocks() const;
    std::string pickle() const;
    static std::unique_ptr<PostfixNotation> unpickle(const std::string &);
};

class push_reaction_exception : public std::runtime_error {
  public:
    explicit push_reaction_exception(const std::string &message)
        : std::runtime_error(message) {}
};

class Synthesis {
    PostfixNotation postfix_notation;
    std::vector<MolSet> stack;

  public:
    Synthesis() = default;
    Synthesis(const PostfixNotation &&postfix_notation,
              const std::vector<MolSet> &&stack)
        : postfix_notation(postfix_notation), stack(stack) {}
    Synthesis(const Synthesis &) = default;
    Synthesis(Synthesis &&) = default;
    Synthesis &operator=(const Synthesis &) = default;
    const PostfixNotation &get_postfix_notation() const;
    const std::vector<MolSet> &get_stack() const;
    void push(const Mol_sptr &);
    void push(const Reaction_sptr &, size_t max_products = 8);
    void push(const Synthesis &);
    void push(const std::shared_ptr<Synthesis> &);
    void push(const std::variant<Mol_sptr, std::shared_ptr<Synthesis>> &);
    MolSet top(std::vector<MolSet>::size_type) const;
    MolSet top() const;
    std::vector<MolSet>::size_type stack_size() const;
    bool is_empty() const;
    size_t count_reactions() const;
    size_t count_building_blocks() const;
    std::string pickle() const;
    static std::unique_ptr<Synthesis> unpickle(const std::string &);
};

typedef std::shared_ptr<Synthesis> Synthesis_sptr;
typedef std::vector<Synthesis_sptr> SynthesisVector;
void save_synthesis_vector(const SynthesisVector &,
                           const std::filesystem::path &);
SynthesisVector *load_synthesis_vector(const std::filesystem::path &path);
} // namespace prexsyn_engine
