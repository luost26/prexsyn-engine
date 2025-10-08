#pragma once

#include <functional>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#include "../chemspace.hpp"
#include "../featurizer/featurizer.hpp"
#include "../utils/assert.hpp"
#include "../utils/logging.hpp"
#include "data_buffer.hpp"

namespace prexsyn_engine {

template <size_t capacity> class DataPipelineV2 {
    size_t num_requested_threads;

    std::shared_ptr<ChemicalSpaceDefinition> csd;
    SynthesisGeneratorOption gen_option;
    using Builder = typename DataBuffer<capacity>::WriteTransaction;
    std::shared_ptr<FeaturizerSet> featurizer;
    std::mt19937::result_type base_seed;

    std::vector<std::jthread> thread_pool;
    DataBuffer<capacity> buffer;

    void worker_fn(int worker_id) {
        auto seed = base_seed + worker_id;
        logger()->info("Thread {} started, seed={}", worker_id, seed);
        SynthesisGenerator generator(csd, gen_option, seed);
        auto st = thread_pool.at(worker_id).get_stop_token();
        while (!st.stop_requested()) {
            auto synthesis = generator.next();
            {
                auto write_txn = buffer.begin_buffered_write();
                (*featurizer)(*synthesis, *write_txn);
                write_txn->commit();
            }
        }
    }

  public:
    DataPipelineV2(size_t num_threads,
                   const std::shared_ptr<ChemicalSpaceDefinition> &csd,
                   const SynthesisGeneratorOption &gen_option,
                   const std::shared_ptr<FeaturizerSet> &featurizer,
                   std::mt19937::result_type base_seed = std::random_device{}())
        : num_requested_threads(num_threads), csd(csd), gen_option(gen_option),
          featurizer(featurizer), base_seed(base_seed) {
        Ensures(num_requested_threads > 0);
    }

    void start() {
        for (size_t i = 0; i < num_requested_threads; ++i) {
            thread_pool.emplace_back(&DataPipelineV2<capacity>::worker_fn, this,
                                     i);
        }
    }

    using GetCallback =
        std::function<void(const std::vector<typename DataBuffer<
                               capacity>::ReadTransaction::ReadEntry> &)>;

    void get(const GetCallback &callback, size_t n = 1) {
        Ensures(!thread_pool.empty());
        auto read_txn = buffer.begin_read(n); // scope lock acquired
        auto entries = read_txn->read_all();
        callback(entries);
    }

    void stop() {
        for (auto &thread : thread_pool) {
            thread.request_stop();
        }
        buffer.clear();
        for (auto &thread : thread_pool) {
            thread.join();
            logger()->info("Thread joined");
        }
        buffer.clear();
    }

    ~DataPipelineV2() { stop(); }
};

} // namespace prexsyn_engine
