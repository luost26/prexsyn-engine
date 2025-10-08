#include "data_buffer.hpp"

#include <gtest/gtest.h>

#include "../utils/logging.hpp"

using namespace prexsyn_engine;

TEST(data_buffer, case1) {
    DataBuffer<16> buffer{};

    for (int i = 100; i < 10000; i += 100) {
        {
            auto write_txn = buffer.begin_buffered_write();
            logger()->info("Starting write transaction. cursor={}",
                           write_txn->cursor());
            write_txn->add("test1", (long)i + 42);
            write_txn->add("test2", std::vector<long>{i + 1, i + 2, i + 3});
            write_txn->commit();
        }

        {
            auto write_txn = buffer.begin_buffered_write();
            logger()->info("Starting write transaction. cursor={}",
                           write_txn->cursor());
            write_txn->add("test1", (long)i + 43);
            write_txn->add("test2", std::vector<long>{i + 4, i + 5, i + 6});
            write_txn->commit();
        }

        {
            auto read_txn = buffer.begin_read(
                2); // 16 is multiple of 2 so no span2 is needed
            logger()->info("Starting read transaction. cursor={}",
                           read_txn->cursor());
            auto entries = read_txn->read_all();

            ASSERT_EQ(entries.size(), 2);
            EXPECT_EQ(entries[0].name, "test1");
            EXPECT_EQ(entries[0].shape, std::vector<long>{});
            EXPECT_EQ(entries[0].dtype, DType::Long);
            EXPECT_EQ(*reinterpret_cast<long *>(entries[0].span1.data()),
                      i + 42);
            EXPECT_EQ(*(reinterpret_cast<long *>(entries[0].span1.data()) + 1),
                      i + 43);

            EXPECT_EQ(entries[1].name, "test2");
            EXPECT_EQ(entries[1].shape, std::vector<long>{3});
            EXPECT_EQ(entries[1].dtype, DType::Long);
            auto vec = reinterpret_cast<long *>(entries[1].span1.data());
            EXPECT_EQ(vec[0], i + 1);
            EXPECT_EQ(vec[1], i + 2);
            EXPECT_EQ(vec[2], i + 3);
            EXPECT_EQ(vec[3], i + 4);
            EXPECT_EQ(vec[4], i + 5);
            EXPECT_EQ(vec[5], i + 6);
        }
    }
}
