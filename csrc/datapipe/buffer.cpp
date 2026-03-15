#include "buffer.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace prexsyn::datapipe {

template <size_t capacity>
    requires(capacity > 0)
DataBuffer<capacity>::DataBuffer(const Schema &schema) : schema_(schema) {
    columns_.reserve(schema.size());
    for (size_t col_index = 0; col_index < schema.size(); ++col_index) {
        columns_.emplace_back(schema[col_index]);
        column_name_to_index_.emplace(schema[col_index].name(), col_index);
    }
}

template <size_t capacity>
    requires(capacity > 0)
std::unique_ptr<WriteRow<capacity>> DataBuffer<capacity>::new_write_row() {
    return std::unique_ptr<WriteRow<capacity>>(new WriteRow<capacity>(*this));
}

template <size_t capacity>
    requires(capacity > 0)
void DataBuffer<capacity>::put(std::unique_ptr<WriteRow<capacity>> row) {
    empty_sem.acquire();
    {
        std::scoped_lock lock(mutex);
        for (size_t col_index = 0; col_index < columns_.size(); ++col_index) {
            auto &column = columns_[col_index];
            auto span = column.raw_span(write_cursor);
            std::copy(row->row_data_[col_index].begin(), row->row_data_[col_index].end(),
                      span.begin());
        }
        write_cursor = (write_cursor + 1) % capacity;
    }
    full_sem.release();
}

template <size_t capacity>
    requires(capacity > 0)
void DataBuffer<capacity>::get(const ReadBatch &batch) {
    if (batch.batch_size > capacity) {
        throw std::runtime_error("Batch size cannot be greater than buffer capacity");
    }
    if (batch.destinations.size() != columns_.size()) {
        throw std::runtime_error("Batch destinations size must match number of columns");
    }

    for (size_t i = 0; i < batch.batch_size; ++i) {
        full_sem.acquire();
    }
    {
        std::scoped_lock lock(mutex);

        size_t start = read_cursor;
        size_t end = read_cursor + batch.batch_size;
        bool use2 = false;
        size_t start1{}, end1{}, start2{}, end2{};
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

        for (size_t col_index = 0; col_index < columns_.size(); ++col_index) {
            auto &column = columns_[col_index];
            auto dst_span = batch.destinations.at(col_index);

            auto src_span1 = column.raw_span(start1, end1);
            auto dst_span1 = dst_span.subspan(0, src_span1.size());
            std::copy(src_span1.begin(), src_span1.end(), dst_span1.begin());

            if (use2) {
                auto src_span2 = column.raw_span(start2, end2);
                auto dst_span2 = dst_span.subspan(src_span1.size(), src_span2.size());
                std::copy(src_span2.begin(), src_span2.end(), dst_span2.begin());
            }
        }

        read_cursor = (read_cursor + batch.batch_size) % capacity;
    }
    for (size_t i = 0; i < batch.batch_size; ++i) {
        empty_sem.release();
    }
}

template <size_t capacity>
    requires(capacity > 0)
void DataBuffer<capacity>::get(const NamedReadBatch &batch) {
    ReadBatch index_batch;
    index_batch.batch_size = batch.batch_size;
    index_batch.destinations.reserve(columns_.size());
    for (size_t col_index = 0; col_index < columns_.size(); ++col_index) {
        const auto &col_name = columns_[col_index].name();
        if (batch.destinations.find(col_name) == batch.destinations.end()) {
            throw std::runtime_error("Batch destinations missing column: " + col_name);
        }
        index_batch.destinations.push_back(batch.destinations.at(col_name));
    }
    get(index_batch);
}

template class DataBuffer<4>;
template class DataBuffer<8192>;
template class DataBuffer<65536>;

} // namespace prexsyn::datapipe
