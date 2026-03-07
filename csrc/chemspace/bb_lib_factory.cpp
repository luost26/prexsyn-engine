#include "bb_lib_factory.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>

#include <GraphMol/FileParsers/MolSupplier.h>

#include "../chemistry/chemistry.hpp"
#include "../utility/logging.hpp"
#include "bb_lib.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<BuildingBlockLibrary>
bb_lib_from_sdf(const std::filesystem::path &path, const BuildingBlockPreprocessor &preprocessor) {
    auto bb_lib = std::make_unique<BuildingBlockLibrary>();

    RDKit::SDMolSupplier supplier(path.string(), true, false, true);
    logger()->info("Starting to load building blocks from SDF: {}", path.string());
    size_t count = 0;
    while (!supplier.atEnd()) {
        try {
            RDKit::ROMOL_SPTR rdkit_mol{supplier.next()};
            auto mol = preprocessor(std::make_shared<Molecule>(rdkit_mol));
            if (!rdkit_mol->hasProp("id")) {
                logger()->warn("'id' property missing in SDF entry: {}", mol->smiles());
                continue;
            }
            bb_lib->add({
                .molecule = std::move(mol),
                .identifier = rdkit_mol->getProp<std::string>("id"),
                .classifications = {},
            });
            count++;
            if (count % 10000 == 0) {
                logger()->info("Loaded {} building blocks ...", count);
            }
        } catch (const MoleculeError &e) {
            logger()->warn("Failed to process SDF entry: {}", e.what());
            continue;
        }
    }
    logger()->info("Done. Total loaded: {}", count);
    return bb_lib;
}

} // namespace prexsyn::chemspace
