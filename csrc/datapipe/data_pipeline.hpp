#pragma once

#include <cstddef>
#include <map>
#include <memory>
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
};

class Worker {
private:
    friend class DataPipeline;

    const DataPipeline &owner_;
    Generator generator_;
    std::jthread thread_;

    Worker(const DataPipeline &owner, size_t seed)
        : owner_(owner), generator_(owner_.chemical_space_, owner_.generator_config_, seed),
          thread_(&Worker::run, this) {
        owner_.logger_->info("Worker thread started with seed {}", seed);
    }

    void request_stop() { thread_.request_stop(); }
    void join() {
        if (thread_.joinable()) {
            thread_.join();
        }
    }

    void run();
};

} // namespace prexsyn::datapipe
