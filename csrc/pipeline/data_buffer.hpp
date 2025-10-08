#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <semaphore>
#include <span>
#include <vector>

#include "../featurizer/builder.hpp"
#include "../utils/assert.hpp"
#include "../utils/logging.hpp"

namespace prexsyn_engine {

template <size_t capacity>
    requires(capacity > 0)
class DataBuffer {
  public:
    class WriteTransaction : public FeatureBuilder {
        DataBuffer &buffer;
        std::map<std::string, std::vector<long>> shape;
        std::map<std::string, DType::Type> dtype;
        std::map<std::string, std::unique_ptr<std::byte[]>> data;
        bool committed = false;

        template <supported_dtype T>
        void init(const std::string &name, const std::vector<long> &shape) {
            if (this->shape.find(name) != this->shape.end()) {
                throw std::runtime_error("Duplicate write to " + name);
            }
            this->shape[name] = shape;
            dtype[name] = DType::get_type<T>();
            size_t numel = 1;
            for (auto dim : shape) {
                numel *= dim;
            }
            data[name] = std::make_unique<std::byte[]>(numel * sizeof(T));
        }

        size_t get_numel(const std::string &name) {
            Ensures(shape.find(name) != shape.end());
            size_t numel = 1;
            for (const auto dim : shape[name]) {
                numel *= dim;
            }
            return numel;
        }

        template <typename T> T *get_pointer(const std::string &name) {
            Ensures(dtype[name] == DType::get_type<T>());
            Ensures(data.find(name) != data.end());
            return reinterpret_cast<T *>(data[name].get());
        }

      public:
        WriteTransaction(DataBuffer &buffer) : buffer(buffer) {}
        WriteTransaction(const WriteTransaction &) = delete;
        WriteTransaction &operator=(const WriteTransaction &) = delete;
        WriteTransaction(WriteTransaction &&) = delete;
        WriteTransaction &operator=(WriteTransaction &&) = delete;

        size_t cursor() const { return buffer.write_cursor; }
        void commit() {
            Ensures(!committed);
            committed = true;
            buffer.empty_sem.acquire();
            {
                std::scoped_lock lock(buffer.mutex);
                for (const auto &[name, _] : shape) {
                    if (!buffer.has(name)) {
                        buffer.init(name, shape[name], dtype[name]);
                    } else {
                        Ensures(buffer.shape[name] == shape[name]);
                        Ensures(buffer.dtype[name] == dtype[name]);
                    }

                    auto item_size =
                        get_numel(name) * DType::get_size(dtype[name]);
                    auto offset = buffer.write_cursor * item_size;
                    auto ptr = buffer.get_raw_pointer(name) + offset;
                    std::memcpy(ptr, data[name].get(), item_size);
                }

                buffer.write_cursor = (buffer.write_cursor + 1) % capacity;
            }
            buffer.full_sem.release();
        }

#define DEFINE_ADD_METHODS(T)                                                  \
    void add(const std::string &name, const T &scalar) {                       \
        std::vector<long> shape = shape_of(scalar);                            \
        init<T>(name, shape);                                                  \
                                                                               \
        get_pointer<T>(name)[0] = scalar;                                      \
    }                                                                          \
                                                                               \
    void add(const std::string &name, const std::vector<T> &vec) {             \
        std::vector<long> shape = shape_of(vec);                               \
        init<T>(name, shape);                                                  \
                                                                               \
        auto numel = get_numel(name);                                          \
        auto ptr = get_pointer<T>(name);                                       \
        for (size_t i = 0; i < numel; ++i) {                                   \
            ptr[i] = vec[i];                                                   \
        }                                                                      \
    }                                                                          \
                                                                               \
    void add(const std::string &name,                                          \
             const std::vector<std::vector<T>> &vec) {                         \
        std::vector<long> shape = shape_of(vec);                               \
        init<T>(name, shape);                                                  \
                                                                               \
        auto n_rows = shape[0];                                                \
        auto n_cols = shape[1];                                                \
        auto ptr = get_pointer<T>(name);                                       \
        for (auto i = 0; i < n_rows; ++i) {                                    \
            for (auto j = 0; j < n_cols; ++j) {                                \
                ptr[i * n_cols + j] = vec[i][j];                               \
            }                                                                  \
        }                                                                      \
    }
        DEFINE_ADD_METHODS(float)
        DEFINE_ADD_METHODS(long)
        DEFINE_ADD_METHODS(bool)
#undef DEFINE_ADD_METHODS
    };

    class ReadTransaction {
        DataBuffer &buffer;
        size_t n;

      public:
        struct ReadEntry {
            std::string name;
            std::vector<long> shape;
            DType::Type dtype;
            std::span<std::byte> span1;
            std::optional<std::span<std::byte>> span2;

            ReadEntry(const std::string &name, const std::vector<long> &shape,
                      DType::Type dtype, std::span<std::byte> span1,
                      std::optional<std::span<std::byte>> span2)
                : name(name), shape(shape), dtype(dtype), span1(span1),
                  span2(span2) {}
        };

        ReadTransaction(DataBuffer &buffer, size_t n = 1)
            : buffer(buffer), n(n) {
            Ensures(n > 0 && n <= capacity);
            for (size_t i = 0; i < n; ++i) {
                buffer.full_sem.acquire();
            }
            buffer.mutex.lock();
        }
        ~ReadTransaction() {
            buffer.read_cursor = (buffer.read_cursor + n) % capacity;
            buffer.mutex.unlock();
            for (size_t i = 0; i < n; ++i) {
                buffer.empty_sem.release();
            }
        }
        ReadTransaction(const ReadTransaction &) = delete;
        ReadTransaction &operator=(const ReadTransaction &) = delete;
        ReadTransaction(ReadTransaction &&) = delete;
        ReadTransaction &operator=(ReadTransaction &&) = delete;

        size_t cursor() const { return buffer.read_cursor; }

        std::vector<ReadEntry> read_all() {
            size_t start = buffer.read_cursor, end = buffer.read_cursor + n;
            bool use2 = false;
            size_t start1, end1, start2, end2;
            if (end >= capacity) {
                start1 = start;
                end1 = capacity;
                use2 = true;
                start2 = 0;
                end2 = end - capacity;
            } else {
                start1 = start;
                end1 = end;
                use2 = false;
                start2 = 0;
                end2 = 0;
            }
            std::vector<ReadEntry> entries;
            for (const auto &[name, shape] : buffer.shape) {
                auto dtype = buffer.dtype[name];
                auto numel = buffer.get_numel(name);
                std::byte *raw_ptr = buffer.get_raw_pointer(name);

                auto offset1 = start1 * numel * DType::get_size(dtype);
                auto size1 = numel * (end1 - start1) * DType::get_size(dtype);
                auto span1 = std::span<std::byte>(raw_ptr + offset1, size1);

                std::optional<std::span<std::byte>> span2 = std::nullopt;
                if (use2) {
                    auto offset2 = start2 * numel * DType::get_size(dtype);
                    auto size2 =
                        numel * (end2 - start2) * DType::get_size(dtype);
                    span2 = std::span<std::byte>(raw_ptr + offset2, size2);
                }
                entries.emplace_back(name, shape, dtype, span1, span2);
            }
            return entries;
        }
    };

  private:
    std::map<std::string, std::vector<long>> shape;
    std::map<std::string, DType::Type> dtype;
    std::map<std::string, std::unique_ptr<std::byte[]>> data;
    std::counting_semaphore<capacity> empty_sem{capacity};
    std::counting_semaphore<capacity> full_sem{0};
    std::mutex mutex;
    size_t write_cursor = 0;
    size_t read_cursor = 0;

    bool has(const std::string &name) { return data.find(name) != data.end(); }

    size_t get_numel(const std::string &name) {
        Ensures(has(name));
        size_t numel = 1;
        for (auto dim : shape[name]) {
            numel *= dim;
        }
        return numel;
    }

    void init(const std::string &name, const std::vector<long> &shape,
              const DType::Type &dtype) {
        Ensures(!has(name));
        logger()->info("Initializing buffer for {}", name);

        this->shape[name] = shape;
        this->dtype[name] = dtype;
        size_t numel = 1;
        for (auto dim : shape) {
            numel *= dim;
        }
        Ensures(numel > 0);
        data[name] = std::make_unique<std::byte[]>(capacity * numel *
                                                   DType::get_size(dtype));
    }

    std::byte *get_raw_pointer(const std::string &name) {
        Ensures(has(name));
        return data[name].get();
    }

    template <typename T> T *get_pointer(const std::string &name) {
        Ensures(DType::get_type<T>() == dtype[name]);
        return reinterpret_cast<T *>(data[name].get());
    }

  public:
    std::unique_ptr<WriteTransaction> begin_buffered_write() {
        return std::make_unique<WriteTransaction>(*this);
    }

    std::unique_ptr<ReadTransaction> begin_read(size_t n = 1) {
        return std::make_unique<ReadTransaction>(*this, n);
    }

    void clear() {
        std::scoped_lock lock(mutex);
        while (full_sem.try_acquire()) {
            empty_sem.release();
        }
        read_cursor = 0;
        write_cursor = 0;
        shape.clear();
        dtype.clear();
        data.clear();
    }
};
} // namespace prexsyn_engine
