#include "chemspace.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <omp.h>

#include "synthesis.hpp"
#include "utils/algorithm.hpp"
#include "utils/assert.hpp"
#include "utils/logging.hpp"

namespace synthesis_backend {

const std::shared_ptr<BuildingBlockList>
ChemicalSpaceDefinition::get_primary_building_blocks() const {
    return primary_building_blocks;
}

const std::shared_ptr<SynthesisVector>
ChemicalSpaceDefinition::get_secondary_building_blocks() const {
    return secondary_building_blocks;
}

const std::shared_ptr<ReactionList>
ChemicalSpaceDefinition::get_reactions() const {
    return reactions;
}

void ChemicalSpaceDefinition::save(const std::filesystem::path &dirpath) const {
    if (!std::filesystem::exists(dirpath)) {
        std::filesystem::create_directories(dirpath);
    }
    primary_building_blocks->save(dirpath / "primary_building_blocks");
    save_synthesis_vector(*secondary_building_blocks,
                          dirpath / "secondary_building_blocks");
    reactions->save(dirpath / "reactions");
    primary_index->save(dirpath / "primary_index");
    secondary_index->save(dirpath / "secondary_index");
}

std::variant<Mol_sptr, Synthesis_sptr>
ChemicalSpaceDefinition::random_building_block(std::mt19937 &rng) const {
    return random_choice(primary_building_blocks->get_building_blocks(),
                         *secondary_building_blocks, rng);
}

std::variant<Mol_sptr, Synthesis_sptr>
ChemicalSpaceDefinition::random_building_block(
    const ReactionReactantIndexTuple &rr_idx, std::mt19937 &rng) const {
    const std::vector<MolecularIndex> &primary_indices =
        primary_index->get_molecular_indices(rr_idx.reaction, rr_idx.reactant);
    const std::vector<MolecularIndex> &secondary_indices =
        secondary_index->get_molecular_indices(rr_idx.reaction,
                                               rr_idx.reactant);

    try {
        bool is_secondary;
        MolecularIndex choice = random_choice(
            primary_indices, secondary_indices, rng, &is_secondary);
        if (is_secondary) {
            return secondary_building_blocks->at(choice);
        } else {
            return primary_building_blocks->get(choice);
        }
    } catch (const std::out_of_range &e) {
        throw no_available_building_blocks(
            "No available building blocks for the given reaction and "
            "reactant index.");
    }
}

std::vector<ReactionReactantIndexTuple>
ChemicalSpaceDefinition::available_reactions(const Mol_sptr &mol) const {
    std::vector<ReactionReactantIndexTuple> result;
    for (ReactionIndex i = 0; i < reactions->get_reactions().size(); ++i) {
        const auto &reaction = reactions->get(i);
        std::vector<ReactantIndex> reactant_indices =
            get_suitable_reactant_indices(reaction, mol);
        for (const auto &reactant_index : reactant_indices) {
            result.push_back({i, reactant_index});
        }
    }
    return result;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::building_blocks_from_sdf(
    const std::filesystem::path &path,
    const BuildingBlockPreprocessingOption &option) {
    this->primary_building_blocks.reset(
        BuildingBlockList::from_sdf(path, option));
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::building_blocks_from_cache(
    const std::filesystem::path &path) {
    this->primary_building_blocks.reset(BuildingBlockList::load(path));
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::reactions_from_txt(
    const std::filesystem::path &path) {
    this->reactions.reset(ReactionList::from_txt(path));
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::reactions_from_cache(
    const std::filesystem::path &path) {
    this->reactions.reset(ReactionList::load(path));
    return *this;
}

ChemicalSpaceDefinitionBuilder &ChemicalSpaceDefinitionBuilder::
    secondary_building_blocks_from_single_reaction() {
    Ensures(reactions != nullptr);
    Ensures(primary_building_blocks != nullptr);

    logger()->info(
        "Generating secondary building blocks from single reactions...");

    secondary_building_blocks.reset(new SynthesisVector());
#pragma omp parallel for
    for (const auto &building_block :
         primary_building_blocks->get_building_blocks()) {
        for (const auto &reaction : reactions->get_reactions()) {
            Synthesis_sptr synthesis(new Synthesis());
            if (reaction->getNumReactantTemplates() != 1) {
                continue;
            }
            try {
                synthesis->push(building_block);
                synthesis->push(reaction);
#pragma omp critical
                {
                    secondary_building_blocks->push_back(synthesis);
                }
            } catch (const push_reaction_exception &e) {
                continue;
            }
        }
    }
    logger()->info("Generated {} secondary building blocks.",
                   secondary_building_blocks->size());

    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::secondary_building_blocks_from_cache(
    const std::filesystem::path &path) {
    this->secondary_building_blocks.reset(load_synthesis_vector(path));
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::no_secondary_building_blocks() {
    this->secondary_building_blocks.reset(new SynthesisVector());
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::build_primary_index() {
    Ensures(primary_building_blocks != nullptr);
    Ensures(reactions != nullptr);

    this->primary_index.reset(new ReactionToMolecular<Mol_sptr>(
        primary_building_blocks->get_building_blocks(),
        reactions->get_reactions()));
    return *this;
}

ChemicalSpaceDefinitionBuilder &
ChemicalSpaceDefinitionBuilder::build_secondary_index() {
    Ensures(secondary_building_blocks != nullptr);
    Ensures(reactions != nullptr);

    this->secondary_index.reset(new ReactionToMolecular<Synthesis_sptr>(
        *secondary_building_blocks, reactions->get_reactions()));
    return *this;
}

ChemicalSpaceDefinitionBuilder &ChemicalSpaceDefinitionBuilder::all_from_cache(
    const std::filesystem::path &dirpath) {

    this->primary_building_blocks.reset(
        BuildingBlockList::load(dirpath / "primary_building_blocks"));
    this->secondary_building_blocks.reset(
        load_synthesis_vector(dirpath / "secondary_building_blocks"));
    this->reactions.reset(ReactionList::load(dirpath / "reactions"));
    this->primary_index.reset(
        ReactionToMolecular<Mol_sptr>::load(dirpath / "primary_index"));
    this->secondary_index.reset(
        ReactionToMolecular<Synthesis_sptr>::load(dirpath / "secondary_index"));
    return *this;
}

std::shared_ptr<ChemicalSpaceDefinition>
ChemicalSpaceDefinitionBuilder::build() {
    Ensures(primary_building_blocks != nullptr);
    Ensures(secondary_building_blocks != nullptr);
    Ensures(reactions != nullptr);
    Ensures(primary_index != nullptr);
    Ensures(secondary_index != nullptr);

    auto csd = std::make_shared<ChemicalSpaceDefinition>();
    csd->primary_building_blocks = std::move(primary_building_blocks);
    csd->secondary_building_blocks = std::move(secondary_building_blocks);
    csd->reactions = std::move(reactions);
    csd->primary_index = std::move(primary_index);
    csd->secondary_index = std::move(secondary_index);

    primary_building_blocks = nullptr;
    secondary_building_blocks = nullptr;
    reactions = nullptr;
    primary_index = nullptr;
    secondary_index = nullptr;
    return csd;
}

SynthesisVector::size_type
ChemicalSpaceDefinitionBuilder::count_secondary_building_blocks() const {
    if (secondary_building_blocks) {
        return secondary_building_blocks->size();
    } else {
        throw std::runtime_error(
            "Secondary building blocks are not initialized.");
    }
}

Synthesis_sptr ChemicalSpaceDefinitionBuilder::get_secondary_building_block(
    SynthesisVector::size_type index) const {
    if (secondary_building_blocks) {
        return secondary_building_blocks->at(index);
    } else {
        throw std::runtime_error(
            "Secondary building blocks are not initialized.");
    }
}

SynthesisGenerator::SynthesisGenerator(
    const std::shared_ptr<ChemicalSpaceDefinition> &csd,
    const SynthesisGeneratorOption &option, std::mt19937::result_type seed)
    : csd(csd), synthesis(new Synthesis()), rng(seed), option(option) {}

SynthesisGenerator::SynthesisGenerator(const SynthesisGenerator &other)
    : csd(other.csd), synthesis(new Synthesis(*other.synthesis)) {}

void SynthesisGenerator::reset_synthesis() { synthesis.reset(new Synthesis()); }

void SynthesisGenerator::init_synthesis() {
    Ensures(synthesis->stack_size() == 0);
    synthesis->push(csd->random_building_block(rng));
}

size_t count_product_atoms(const Synthesis &synthesis) {
    size_t count = 0;
    for (const auto &mol : synthesis.top()) {
        count = std::max(count, (size_t)mol->getNumAtoms());
    }
    return count;
}

bool SynthesisGenerator::needs_reinit() const {
    if (synthesis->count_reactions() >= option.num_reactions_cutoff) {
        return true;
    } else if (count_product_atoms(*synthesis) >=
               option.num_product_atoms_cutoff) {
        return true;
    }
    return false;
}

Synthesis_sptr SynthesisGenerator::next() {
    if (synthesis->stack_size() == 0) {
        init_synthesis();
    } else {
        try {
            auto available =
                csd->available_reactions(random_choice(synthesis->top(), rng));
            if (available.empty()) {
                throw no_available_reactions(
                    "No available reactions for the given synthesis.");
            }
            auto rr_idx = random_choice(available, rng);
            auto total_reactants =
                csd->reactions->get(rr_idx.reaction)->getNumReactantTemplates();
            for (size_t other_reactant_idx = 0;
                 other_reactant_idx < total_reactants; ++other_reactant_idx) {
                if (other_reactant_idx != rr_idx.reactant) {
                    auto reactant = csd->random_building_block(
                        {rr_idx.reaction, other_reactant_idx}, rng);
                    synthesis->push(reactant);
                }
            }
            synthesis->push(csd->reactions->get(rr_idx.reaction));
        } catch (const push_reaction_exception &e) {
            logger()->debug("Push reaction exception: {}", e.what());
            reset_synthesis();
            init_synthesis();
        } catch (const no_available_building_blocks &e) {
            logger()->debug("No available building blocks: {}", e.what());
            reset_synthesis();
            init_synthesis();
        } catch (const no_available_reactions &e) {
            logger()->debug("No available reactions: {}", e.what());
            reset_synthesis();
            init_synthesis();
        }
    }
    auto out = std::make_shared<Synthesis>(*synthesis);
    if (needs_reinit()) {
        reset_synthesis();
    }
    Ensures(out->stack_size() == 1);
    return out;
}
} // namespace synthesis_backend
