#include "data_pipeline.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <stop_token>
#include <string>
#include <utility>
#include <vector>

#include "../utility/logging.hpp"
#include "buffer.hpp"
#include "generator.hpp"

namespace prexsyn::datapipe {

DataPipeline::DataPipeline(const std::shared_ptr<ChemicalSpace> &cs,
                           const std::map<std::string, std::shared_ptr<MoleculeDescriptor>> &md,
                           const std::map<std::string, std::shared_ptr<SynthesisDescriptor>> &sd,
                           const GeneratorConfig &generator_config)
    : chemical_space_(cs), generator_config_(generator_config), molecule_descriptors_(md),
      synthesis_descriptors_(sd), logger_(create_logger("DataPipeline")) {
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

void DataPipeline::start_workers(const std::vector<size_t> &seeds) {
    for (const auto &seed : seeds) {
        std::unique_ptr<Worker> worker{new Worker(*this, seed)};
        workers_.push_back(std::move(worker));
    }
}

void DataPipeline::stop_workers() {
    for (auto &worker : workers_) {
        worker->request_stop();
    }
    for (auto &worker : workers_) {
        worker->join();
    }
    workers_.clear();
}

void DataPipeline::get(const NamedReadBatch &batch) { buffer_->get(batch); }

void Worker::run() {
    while (!thread_.get_stop_token().stop_requested()) {
        auto [synthesis, product] = generator_.next_with_product();
        auto data_row = owner_.buffer_->new_write_row();
        for (const auto &[name, fn] : owner_.synthesis_descriptors_) {
            auto dest_span = data_row->data(name);
            (*fn)(*synthesis, dest_span);
        }
        for (const auto &[name, fn] : owner_.molecule_descriptors_) {
            auto dest_span = data_row->data(name);
            (*fn)(*product, dest_span);
        }
        owner_.buffer_->put(std::move(data_row));
    }
}

} // namespace prexsyn::datapipe
