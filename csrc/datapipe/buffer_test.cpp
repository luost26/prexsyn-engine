#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "buffer.hpp"

namespace {

using prexsyn::DataType;
using prexsyn::datapipe::ColumnDef;
using prexsyn::datapipe::DataBuffer;
using prexsyn::datapipe::ReadBatch;

constexpr size_t kTestCapacity = 4;

std::vector<ColumnDef> make_schema() {
    return {
        ColumnDef("scores", {2}, DataType::T::float32),
        ColumnDef("ids", {1}, DataType::T::int64),
    };
}

template <typename T> std::span<std::byte> as_writable_bytes(std::vector<T> &values) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    return {reinterpret_cast<std::byte *>(values.data()), values.size() * sizeof(T)};
}

void write_row(DataBuffer<kTestCapacity> &buffer, std::array<float, 2> scores, std::int64_t id) {
    auto row = buffer.new_write_row();
    auto score_span = row->data<float>("scores");
    auto id_span = row->data<std::int64_t>("ids");

    std::copy(scores.begin(), scores.end(), score_span.begin());
    id_span[0] = id;

    buffer.put(std::move(row));
}

void expect_row(const std::vector<float> &scores, const std::vector<std::int64_t> &ids,
                size_t row_index, std::array<float, 2> expected_scores, std::int64_t expected_id) {
    EXPECT_FLOAT_EQ(scores.at(row_index * 2), expected_scores[0]);
    EXPECT_FLOAT_EQ(scores.at((row_index * 2) + 1), expected_scores[1]);
    EXPECT_EQ(ids.at(row_index), expected_id);
}

} // namespace

TEST(DataBufferTest, NewWriteRowProvidesTypedColumnViews) {
    DataBuffer<kTestCapacity> buffer(make_schema());

    const auto &column_name_to_index = buffer.column_name_to_index();
    ASSERT_EQ(column_name_to_index.size(), 2);
    EXPECT_EQ(column_name_to_index.at("scores"), 0U);
    EXPECT_EQ(column_name_to_index.at("ids"), 1U);

    auto row = buffer.new_write_row();
    EXPECT_EQ(row->data<float>("scores").size(), 2U);
    EXPECT_EQ(row->data<std::int64_t>("ids").size(), 1U);
    EXPECT_THROW(
        {
            const auto unused = row->data<std::int64_t>("scores");
            (void)unused;
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            const auto unused = row->data<float>("missing");
            (void)unused;
        },
        std::out_of_range);
}

TEST(DataBufferTest, PutAndGetRoundTripMultipleRows) {
    DataBuffer<kTestCapacity> buffer(make_schema());

    write_row(buffer, {1.0F, 1.5F}, 10);
    write_row(buffer, {2.0F, 2.5F}, 20);

    std::vector<float> scores(4);
    std::vector<std::int64_t> ids(2);
    buffer.get(ReadBatch{2, {as_writable_bytes(scores), as_writable_bytes(ids)}});

    expect_row(scores, ids, 0, {1.0F, 1.5F}, 10);
    expect_row(scores, ids, 1, {2.0F, 2.5F}, 20);
}

TEST(DataBufferTest, GetMaintainsOrderAcrossRingWraparound) {
    DataBuffer<kTestCapacity> buffer(make_schema());

    write_row(buffer, {0.0F, 0.5F}, 100);
    write_row(buffer, {1.0F, 1.5F}, 101);
    write_row(buffer, {2.0F, 2.5F}, 102);
    write_row(buffer, {3.0F, 3.5F}, 103);

    std::vector<float> first_scores(6);
    std::vector<std::int64_t> first_ids(3);
    buffer.get(ReadBatch{3, {as_writable_bytes(first_scores), as_writable_bytes(first_ids)}});

    expect_row(first_scores, first_ids, 0, {0.0F, 0.5F}, 100);
    expect_row(first_scores, first_ids, 1, {1.0F, 1.5F}, 101);
    expect_row(first_scores, first_ids, 2, {2.0F, 2.5F}, 102);

    write_row(buffer, {4.0F, 4.5F}, 104);
    write_row(buffer, {5.0F, 5.5F}, 105);
    write_row(buffer, {6.0F, 6.5F}, 106);

    std::vector<float> wrapped_scores(8);
    std::vector<std::int64_t> wrapped_ids(4);
    buffer.get(ReadBatch{4, {as_writable_bytes(wrapped_scores), as_writable_bytes(wrapped_ids)}});

    expect_row(wrapped_scores, wrapped_ids, 0, {3.0F, 3.5F}, 103);
    expect_row(wrapped_scores, wrapped_ids, 1, {4.0F, 4.5F}, 104);
    expect_row(wrapped_scores, wrapped_ids, 2, {5.0F, 5.5F}, 105);
    expect_row(wrapped_scores, wrapped_ids, 3, {6.0F, 6.5F}, 106);
}

TEST(DataBufferTest, GetThrowsWhenBatchExceedsCapacity) {
    DataBuffer<kTestCapacity> buffer(make_schema());

    std::vector<float> scores((kTestCapacity + 1) * 2);
    std::vector<std::int64_t> ids(kTestCapacity + 1);
    EXPECT_THROW(
        {
            buffer.get(
                ReadBatch{kTestCapacity + 1, {as_writable_bytes(scores), as_writable_bytes(ids)}});
        },
        std::runtime_error);
}
