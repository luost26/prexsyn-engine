#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"
#include "../utility/logging.hpp"
#include "buffer.hpp"
#include "generator.hpp"

namespace prexsyn::datapipe {

using chemspace::ChemicalSpace;
using descriptor::MoleculeDescriptor;
using descriptor::SynthesisDescriptor;

class Worker;

class DataPipeline {
private:
    std::shared_ptr<ChemicalSpace> chemical_space_;
    GeneratorConfig generator_config_;

    std::map<std::string, std::shared_ptr<MoleculeDescriptor>> molecule_descriptors_;
    std::map<std::string, std::shared_ptr<SynthesisDescriptor>> synthesis_descriptors_;
    std::unique_ptr<DataBuffer<8192>> buffer_;

    std::shared_ptr<Logger> logger_;

    std::vector<std::unique_ptr<Worker>> workers_;
    friend class Worker;

public:
    DataPipeline(const std::shared_ptr<ChemicalSpace> &,
                 const std::map<std::string, std::shared_ptr<MoleculeDescriptor>> &,
                 const std::map<std::string, std::shared_ptr<SynthesisDescriptor>> &,
                 const GeneratorConfig & = Generator::default_config);

    void start_workers(const std::vector<size_t> &seeds);
    void stop_workers();

    void get(const NamedReadBatch &batch);

    struct Batch {
        const DataPipeline &pipeline;
        size_t batch_size;
        std::map<std::string, std::vector<std::byte>> data;

        template <typename T> std::span<const T> get(const std::string &name) const {
            if (!data.contains(name)) {
                throw std::runtime_error("missing: " + name);
            }
            auto col_idx = pipeline.buffer_->column_name_to_index().at(name);
            auto col_def = pipeline.buffer_->schema().at(col_idx);
            if (DataType::get_dtype<T>() != col_def.dtype()) {
                throw std::runtime_error("type mismatch: " + name);
            }
            const auto &bytes = data.at(name);
            return std::span<const T>(reinterpret_cast<const T *>(bytes.data()),
                                      bytes.size() / sizeof(T));
        }
    };
    Batch get(size_t batch_size);
};

class Worker {
private:
    friend class DataPipeline;

    const DataPipeline &owner_;
    const size_t seed_;
    Generator generator_;
    std::jthread thread_;

    Worker(const DataPipeline &owner, size_t seed);

    void run();
    void request_stop();
    void join();
};

} // namespace prexsyn::datapipe
