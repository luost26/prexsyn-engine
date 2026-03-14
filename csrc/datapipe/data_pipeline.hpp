#pragma once

#include <map>
#include <memory>
#include <string>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"
#include "buffer.hpp"

namespace prexsyn::datapipe {

using chemspace::ChemicalSpace;
using descriptor::MoleculeDescriptor;
using descriptor::SynthesisDescriptor;

class DataPipeline {
private:
    std::shared_ptr<ChemicalSpace> chemical_space_;

    std::map<std::string, std::shared_ptr<MoleculeDescriptor>> molecule_descriptors_;
    std::map<std::string, std::shared_ptr<SynthesisDescriptor>> synthesis_descriptors_;
    std::unique_ptr<DataBuffer<8192>> buffer_;

public:
    DataPipeline(const std::shared_ptr<ChemicalSpace> &,
                 const std::map<std::string, std::shared_ptr<MoleculeDescriptor>> &,
                 const std::map<std::string, std::shared_ptr<SynthesisDescriptor>> &);
};

} // namespace prexsyn::datapipe
