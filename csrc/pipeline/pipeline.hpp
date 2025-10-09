#pragma once

#include <memory>
#include <mutex>
#include <queue>
#include <semaphore>
#include <thread>
#include <tuple>

#include "../chemspace/chemspace.hpp"
#include "../featurizer/featurizer.hpp"
#include "../synthesis.hpp"
#include "../utils/assert.hpp"
#include "../utils/logging.hpp"

namespace prexsyn_engine {

const int MAX_BUFFER_SIZE = 4096;

template <typename Builder> class DataPipeline {
    size_t num_requested_threads;
    std::shared_ptr<ChemicalSpaceDefinition> csd;
    SynthesisGeneratorOption gen_option;
    std::shared_ptr<FeaturizerSet> featurizer;

    std::vector<std::jthread> thread_pool;

    std::queue<std::tuple<Synthesis_sptr, std::shared_ptr<Builder>>> queue;
    std::counting_semaphore<MAX_BUFFER_SIZE> empty_sem{MAX_BUFFER_SIZE};
    std::counting_semaphore<MAX_BUFFER_SIZE> full_sem{0};
    std::mutex buffer_mutex;

    void worker_fn(int worker_id) {
        logger()->info("Thread {} started", worker_id);
        SynthesisGenerator generator(csd, gen_option);
        auto st = thread_pool.at(worker_id).get_stop_token();
        while (!st.stop_requested()) {
            auto synthesis = generator.next();
            auto feature_builder = std::make_shared<Builder>();
            (*featurizer)(*synthesis, feature_builder->erase_type());
            empty_sem.acquire();
            {
                std::scoped_lock lock(buffer_mutex);
                queue.push(std::make_tuple(synthesis, feature_builder));
            }
            full_sem.release();
        }
    }

    void clean_buffer() {
        while (full_sem.try_acquire()) {
            std::scoped_lock lock(buffer_mutex);
            queue.pop();
            empty_sem.release();
        }
    }

  public:
    DataPipeline(size_t num_threads,
                 const std::shared_ptr<ChemicalSpaceDefinition> &csd,
                 const SynthesisGeneratorOption &gen_option,
                 const std::shared_ptr<FeaturizerSet> &featurizer)
        : num_requested_threads(num_threads), csd(csd), gen_option(gen_option),
          featurizer(featurizer) {
        Ensures(num_requested_threads > 0);
    }

    void start() {
        for (size_t i = 0; i < num_requested_threads; ++i) {
            thread_pool.emplace_back(&DataPipeline<Builder>::worker_fn, this,
                                     i);
        }
    }

    std::tuple<Synthesis_sptr, std::shared_ptr<Builder>> get() {
        Ensures(!thread_pool.empty());
        full_sem.acquire();
        std::tuple<Synthesis_sptr, std::shared_ptr<Builder>> result;
        {
            std::scoped_lock lock(buffer_mutex);
            result = queue.front();
            queue.pop();
        }
        empty_sem.release();
        return result;
    }

    void stop() {
        for (auto &thread : thread_pool) {
            thread.request_stop();
        }
        clean_buffer();
        for (auto &thread : thread_pool) {
            thread.join();
            logger()->info("Thread joined");
        }
        clean_buffer();
    }

    ~DataPipeline() { stop(); }
};

} // namespace prexsyn_engine
