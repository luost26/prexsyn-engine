#include "data_pipeline.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <span>
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

    logger_->info("DataPipeline initialized");
    logger_->info("Buffer schema:");
    for (const auto &col_def : buffer_->schema()) {
        logger_->info("- {}: shape {}, dtype {}", col_def.name(), col_def.shape(),
                      DataType::to_string(col_def.dtype()));
    }
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

DataPipeline::Batch DataPipeline::get(size_t batch_size) {
    Batch batch{.pipeline = *this, .batch_size = batch_size, .data = {}};
    for (const auto &col_def : buffer_->schema()) {
        batch.data[col_def.name()] = std::vector<std::byte>(batch_size * col_def.size_in_bytes());
    }

    NamedReadBatch read_batch{.batch_size = batch_size, .destinations = {}};
    for (const auto &col_def : buffer_->schema()) {
        read_batch.add(col_def.name(), std::span<std::byte>(batch.data[col_def.name()]));
    }

    get(read_batch);
    return batch;
}

Worker::Worker(const DataPipeline &owner, size_t seed)
    : owner_(owner), seed_(seed),
      generator_(owner_.chemical_space_, owner_.generator_config_, seed),
      thread_(&Worker::run, this) {
    owner_.logger_->info("Worker[seed={}] started", seed);
}

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

void Worker::request_stop() { thread_.request_stop(); }

void Worker::join() {
    if (thread_.joinable()) {
        thread_.join();
    }
    owner_.logger_->info("Worker[seed={}] joined", seed_);
}

} // namespace prexsyn::datapipe
