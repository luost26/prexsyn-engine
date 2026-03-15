#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <mutex>
#include <semaphore>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

// IWYU pragma: begin_exports
#include "../utility/data_type.hpp"
// IWYU pragma: end_exports

namespace prexsyn::datapipe {

class ColumnDef {
private:
    std::string name_;
    std::vector<size_t> shape_;
    DataType::T dtype_;

public:
    ColumnDef(const std::string &name, const std::vector<size_t> &shape, DataType::T dtype)
        : name_(name), shape_(shape), dtype_(dtype) {}

    const std::string &name() const { return name_; }
    const std::vector<size_t> &shape() const { return shape_; }
    DataType::T dtype() const { return dtype_; }

    size_t num_elements() const {
        size_t num = 1;
        for (size_t dim : shape_) {
            num *= dim;
        }
        return num;
    }

    size_t size_in_bytes() const { return num_elements() * DataType::get_size(dtype_); }
};

template <size_t capacity>
    requires(capacity > 0)
class Column : public ColumnDef {
private:
    std::vector<std::byte> data_;

public:
    Column(const ColumnDef &def) : ColumnDef(def) { data_.resize(def.size_in_bytes() * capacity); }

    std::span<std::byte> raw_span(size_t row_index) {
        size_t byte_offset = row_index * size_in_bytes();
        return {data_.data() + byte_offset, size_in_bytes()};
    }

    std::span<std::byte> raw_span(size_t row_start, size_t row_end) {
        size_t byte_offset = row_start * size_in_bytes();
        size_t byte_size = (row_end - row_start) * size_in_bytes();
        return {data_.data() + byte_offset, byte_size};
    }
};

template <size_t capacity>
    requires(capacity > 0)
class WriteRow;

struct ReadBatch;
struct NamedReadBatch;

template <size_t capacity>
    requires(capacity > 0)
class DataBuffer {
public:
    using Schema = std::vector<ColumnDef>;

private:
    std::vector<Column<capacity>> columns_;
    Schema schema_;
    std::map<std::string, size_t> column_name_to_index_;

    std::counting_semaphore<capacity> empty_sem{capacity};
    std::counting_semaphore<capacity> full_sem{0};
    std::mutex mutex;
    size_t write_cursor = 0;
    size_t read_cursor = 0;

    friend class WriteRow<capacity>;
    friend struct ReadBatch;

public:
    DataBuffer(const std::vector<ColumnDef> &schema);
    const std::vector<ColumnDef> &schema() const { return schema_; }
    const auto &column_name_to_index() const { return column_name_to_index_; }

    std::unique_ptr<WriteRow<capacity>> new_write_row();
    void put(std::unique_ptr<WriteRow<capacity>> row);
    void get(const ReadBatch &batch);
    void get(const NamedReadBatch &batch);
};

template <size_t capacity>
    requires(capacity > 0)
class WriteRow {
private:
    const DataBuffer<capacity> &buffer_;
    std::vector<std::vector<std::byte>> row_data_;

    WriteRow(const DataBuffer<capacity> &buffer) : buffer_(buffer) {
        row_data_.resize(buffer.columns_.size());
        for (size_t i = 0; i < buffer.columns_.size(); ++i) {
            row_data_[i].resize(buffer.columns_[i].size_in_bytes());
        }
    }

    friend class DataBuffer<capacity>;

public:
    std::span<std::byte> data(const std::string &name) {
        auto index = buffer_.column_name_to_index_.at(name);
        return row_data_[index];
    }

    template <typename T> std::span<T> data(const std::string &name) {
        auto index = buffer_.column_name_to_index_.at(name);
        const ColumnDef &def = buffer_.columns_[index];
        if (DataType::get_dtype<T>() != def.dtype()) {
            throw std::runtime_error("Data type mismatch");
        }
        auto &data = row_data_[index];
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return std::span<T>(reinterpret_cast<T *>(data.data()), data.size() / sizeof(T));
    }
};

struct ReadBatch {
    size_t batch_size;
    std::vector<std::span<std::byte>> destinations;

    void add(std::span<std::byte> dest) { destinations.emplace_back(dest); }

    template <SupportedDataType T> void add(std::span<T> destination) {
        destinations.emplace_back(std::as_writable_bytes(destination));
    }
};

struct NamedReadBatch {
    size_t batch_size;
    std::map<std::string, std::span<std::byte>> destinations;

    void add(const std::string &name, std::span<std::byte> dest) { destinations[name] = dest; }

    template <SupportedDataType T> void add(const std::string &name, std::span<T> destination) {
        destinations[name] = std::as_writable_bytes(destination);
    }
};

} // namespace prexsyn::datapipe
