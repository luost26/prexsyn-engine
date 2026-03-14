#include "data_pipeline.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "buffer.hpp"

namespace prexsyn::datapipe {

DataPipeline::DataPipeline(const std::shared_ptr<ChemicalSpace> &cs,
                           const std::map<std::string, std::shared_ptr<MoleculeDescriptor>> &md,
                           const std::map<std::string, std::shared_ptr<SynthesisDescriptor>> &sd)
    : chemical_space_(cs), molecule_descriptors_(md), synthesis_descriptors_(sd) {
    std::vector<ColumnDef> column_defs;
    column_defs.reserve(molecule_descriptors_.size() + synthesis_descriptors_.size());

    for (const auto &[name, descriptor] : molecule_descriptors_) {
        column_defs.emplace_back(name, descriptor->size(), descriptor->dtype());
    }
    for (const auto &[name, descriptor] : synthesis_descriptors_) {
        column_defs.emplace_back(name, descriptor->size(), descriptor->dtype());
    }

    buffer_ = std::make_unique<DataBuffer<8192>>(column_defs);
}

} // namespace prexsyn::datapipe
