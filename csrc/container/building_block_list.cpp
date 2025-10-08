#include "building_block_list.hpp"

#include <stdexcept>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/RWMol.h>
#include <string>

#include "../utils/logging.hpp"

namespace prexsyn_engine {

template <class Archive>
void BuildingBlockPreprocessingOption::serialize(Archive &ar,
                                                 const unsigned int) {
    ar & remove_Hs;
    ar & largest_fragment_only;
}

class molecular_operation_error : public std::runtime_error {
  public:
    explicit molecular_operation_error(const std::string &message)
        : std::runtime_error(message) {}
};

void guarded_molecular_operation(
    Mol_sptr &mol,
    const std::function<
        Mol_sptr::element_type *(const Mol_sptr::element_type &)> &operation,
    const std::string &name = "") {
    if (!mol) {
        throw molecular_operation_error("Null molecule passed to operation (" +
                                        name + ")");
    }
    auto ptr = operation(*mol);
    if (!ptr) {
        throw molecular_operation_error(
            "Operation failed on molecule returning nullptr (" + name + ")");
    }
    mol.reset(ptr);
}

std::optional<Mol_sptr>
BuildingBlockList::preprocess(const Mol_sptr &mol) const {
    Mol_sptr processed_mol = mol;
    try {
        if (option.largest_fragment_only) {
            guarded_molecular_operation(
                processed_mol,
                [&](const RDKit::ROMol &m) {
                    RDKit::MolStandardize::LargestFragmentChooser
                        largest_fragment_chooser;
                    return largest_fragment_chooser.choose(m);
                },
                "largest_fragment_only");
        }
        if (option.remove_Hs) {
            guarded_molecular_operation(
                processed_mol,
                [&](const RDKit::ROMol &m) {
                    return RDKit::MolOps::removeHs(m);
                },
                "remove_Hs");
        }
        return processed_mol;
    } catch (const molecular_operation_error &e) {
        logger()->warn("BuildingBlockList: Preprocessing operation failed: {}",
                       e.what());
        return std::nullopt;
    }
}

std::optional<BuildingBlockIndex>
BuildingBlockList::preprocess_and_add(const Mol_sptr &mol, int original_index) {
    auto processed_mol = preprocess(mol);
    if (processed_mol) {
        auto next_index = building_blocks.size();
        processed_mol.value()->setProp<int>("original_index", original_index);
        processed_mol.value()->setProp<int>("building_block_index", next_index);
        building_blocks.push_back(processed_mol.value());
        return next_index;
    }
    logger()->warn("BuildingBlockList: Failed to preprocess molecule "
                   "(original_index={}): \"{}\", ignoring it",
                   original_index, mol);
    return std::nullopt;
}

BuildingBlockList::BuildingBlockList(
    const MolVector &molecules, const BuildingBlockPreprocessingOption &option)
    : option(option) {
    building_blocks.reserve(molecules.size());

    int original_index = 0;
    for (const auto &mol : molecules) {
        preprocess_and_add(mol, original_index++);
    }
}

BuildingBlockList::BuildingBlockList(
    RDKit::MolSupplier &supplier,
    const BuildingBlockPreprocessingOption &option)
    : option(option) {
    int original_index = 0;
    while (!supplier.atEnd()) {
        Mol_sptr mol(supplier.next());
        logger()->trace("Loading building block #{}: {}", original_index, mol);
        preprocess_and_add(mol, original_index++);

        if (original_index % 10000 == 0) {
            logger()->info(
                "BuildingBlockList: {} building blocks processed, {} loaded",
                original_index, building_blocks.size());
        }
    }
    logger()->info("BuildingBlockList: Finished. {} building blocks loaded",
                   building_blocks.size());
}

BuildingBlockList *
BuildingBlockList::from_sdf(const std::filesystem::path &path,
                            const BuildingBlockPreprocessingOption &option) {
    logger()->info("Loading building blocks from SDF file: {}", path.string());
    RDKit::SDMolSupplier supplier(path.string(), true, false, true);
    return new BuildingBlockList(supplier, option);
}

const MolVector &BuildingBlockList::get_building_blocks() const {
    return building_blocks;
}

void BuildingBlockList::save(const std::filesystem::path &path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open file for writing");
    }
    boost::archive::binary_oarchive oa(ofs);
    oa << option << building_blocks.size();
    for (const auto &mol : building_blocks) {
        std::string pickle;
        RDKit::MolPickler::pickleMol(*mol, pickle, RDKit::PicklerOps::AllProps);
        oa << pickle;
    }
    ofs.close();
}

BuildingBlockList *BuildingBlockList::load(const std::filesystem::path &path) {
    logger()->info("Loading building blocks from cache: {}", path.string());
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading");
    }
    boost::archive::binary_iarchive ia(ifs);
    size_t num_building_blocks;
    BuildingBlockPreprocessingOption option;
    ia >> option;
    BuildingBlockList *object = new BuildingBlockList({}, option);
    ia >> num_building_blocks;
    object->building_blocks.reserve(num_building_blocks);
    for (size_t i = 0; i < num_building_blocks; ++i) {
        std::string pickle;
        ia >> pickle;
        Mol_sptr mol(new RDKit::ROMol());
        RDKit::MolPickler::molFromPickle(pickle, mol.get(),
                                         RDKit::PicklerOps::AllProps);
        object->building_blocks.push_back(mol);
    }
    ifs.close();
    logger()->info("BuildingBlockList: {} building blocks loaded from cache",
                   object->building_blocks.size());
    return object;
}

Mol_sptr BuildingBlockList::get(size_t index) const {
    if (index >= building_blocks.size()) {
        throw std::out_of_range(
            "Index out of range in BuildingBlockList::get " +
            std::to_string(index) +
            " >= " + std::to_string(building_blocks.size()));
    }
    return building_blocks[index];
}

BuildingBlockIndex BuildingBlockList::size() const {
    return building_blocks.size();
}
} // namespace prexsyn_engine
