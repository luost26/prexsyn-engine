#include "bb_lib_factory.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolSupplier.v1API.h>
#include <csv.hpp>

#include "../chemistry/chemistry.hpp"
#include "../utility/logging.hpp"
#include "bb_lib.hpp"

namespace prexsyn::chemspace {

struct identifier_deduplicator {
    std::multiset<std::string> seen;

    std::string operator()(const std::string &identifier) {
        auto count = seen.count(identifier);
        seen.insert(identifier);
        if (count == 0) {
            return identifier;
        } else {
            auto new_id = identifier + "-" + std::to_string(count);
            logger()->warn("Duplicate identifier detected: {}. Renaming to {}", identifier, new_id);
            return new_id;
        }
    }
};

std::unique_ptr<BuildingBlockLibrary>
bb_lib_from_sdf(const std::filesystem::path &path, const BuildingBlockPreprocessor &preprocessor) {
    auto bb_lib = std::make_unique<BuildingBlockLibrary>();
    identifier_deduplicator deduplicator;

    RDKit::SDMolSupplier supplier(path.string(), true, false, true);
    logger()->info("Starting to load building blocks from SDF: {}", path.string());
    size_t count = 0;
    while (!supplier.atEnd()) {
        try {
            RDKit::ROMOL_SPTR rdkit_mol{supplier.next()};
            auto mol = preprocessor(std::make_shared<Molecule>(rdkit_mol));
            if (!rdkit_mol->hasProp("id")) {
                logger()->warn("'id' property missing: {}", mol->smiles());
                continue;
            }
            bb_lib->add({
                .molecule = std::move(mol),
                .identifier = deduplicator(rdkit_mol->getProp<std::string>("id")),
                .labels = {},
            });
            count++;
            if (count % 10000 == 0) {
                logger()->info("Loaded {} building blocks ...", count);
            }
        } catch (const MoleculeError &e) {
            logger()->warn("MoleculeError: {}", e.what());
            continue;
        } catch (const BuildingBlockLibraryError &e) {
            logger()->warn("BuildingBlockLibraryError: {}", e.what());
            continue;
        }
    }
    logger()->info("Done. Loaded: {}", count);
    return bb_lib;
}

std::unique_ptr<BuildingBlockLibrary>
bb_lib_from_csv(const std::filesystem::path &path, const BuildingBlockCSVConfig &config,
                const BuildingBlockPreprocessor &preprocessor) {
    auto bb_lib = std::make_unique<BuildingBlockLibrary>();
    identifier_deduplicator deduplicator;

    csv::CSVReader reader(path.string());
    logger()->info("Starting to load building blocks from CSV: {}", path.string());
    size_t count = 0, rowno = 0;
    for (auto &row : reader) {
        rowno++;
        try {
            std::string identifier, smiles;
            if (!row[config.identifier_column].try_get(identifier) ||
                !row[config.smiles_column].try_get(smiles)) {
                logger()->warn("Missing required columns at row {}", rowno);
                continue;
            }
            auto mol = preprocessor(Molecule::from_smiles(smiles));
            bb_lib->add({
                .molecule = std::move(mol),
                .identifier = deduplicator(identifier),
                .labels = {},
            });
            count++;
            if (count % 10000 == 0) {
                logger()->info("Loaded {} building blocks ...", count);
            }
        } catch (const MoleculeError &e) {
            logger()->warn("MoleculeError: {}", e.what());
            continue;
        } catch (const BuildingBlockLibraryError &e) {
            logger()->warn("BuildingBlockLibraryError: {}", e.what());
            continue;
        }
    }
    logger()->info("Done. Loaded: {}", count);
    return bb_lib;
}

} // namespace prexsyn::chemspace
