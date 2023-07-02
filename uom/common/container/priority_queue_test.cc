// Copyright @2023 UOM project. All rights reserved.

#include "uom/common/container/priority_queue.h"

#include <array>

#include "gtest/gtest.h"

namespace uom {

TEST(PriorityQueue, Basic) {
  PriorityQueue<int> queue;
  queue.reserve(4);
  EXPECT_EQ(0, queue.size());
  queue.emplace(3);
  queue.emplace(1);
  queue.emplace(2);
  queue.emplace(4);
  std::array expected{1, 2, 3, 4};
  for (int i = 0; i < 4; ++i) {
    const int expected_val = expected[i];
    ASSERT_FALSE(queue.empty());
    ASSERT_EQ(4 - i, queue.size());
    ASSERT_EQ(expected_val, queue.top());
    ASSERT_EQ(expected_val, queue.pop());
  }
  ASSERT_TRUE(queue.empty());
}

}  // namespace uom
