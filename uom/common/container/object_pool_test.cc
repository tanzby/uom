// Copyright @2023 UOM project. All rights reserved.

#include "uom/common/container/object_pool.h"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace uom {

struct Foo {
  int a = 0;
  int b = 0;
};

TEST(ObjectPoolTest, Basic) {
  ObjectPool<Foo> pool;

  // Create 2 objects. Then, 2 objects are recycled.
  std::vector<Foo*> expected_addresse;
  {
    auto object1 = pool.Acquire();
    auto object2 = pool.Acquire();
    ASSERT_NE(object1, nullptr);
    ASSERT_NE(object2, nullptr);
    expected_addresse.emplace_back(object1.get());
    expected_addresse.emplace_back(object2.get());
    ASSERT_EQ(0, pool.GetNumCachedObjects());
    ASSERT_EQ(2, pool.GetNumAllocatedObjects());
  }

  ASSERT_EQ(2, pool.GetNumCachedObjects());

  // Acquire again, expect same data address.
  {
    auto object1 = pool.Acquire();
    auto object2 = pool.Acquire();
    const std::vector<Foo*> actual_addresse{object1.get(), object2.get()};
    ASSERT_THAT(actual_addresse, testing::UnorderedElementsAreArray(expected_addresse));
    ASSERT_EQ(0, pool.GetNumCachedObjects());
  }
}

TEST(ObjectPoolTest, Reserve) {
  ObjectPool<Foo> pool;
  pool.Reserve(10);
  ASSERT_EQ(10, pool.GetNumAllocatedObjects());
  ASSERT_EQ(10, pool.GetNumCachedObjects());

  using FooPtr = ObjectPool<Foo>::Ptr;
  std::vector<FooPtr> objects;
  objects.reserve(objects.size());
  for (int i = 0; i < pool.GetNumAllocatedObjects(); ++i) {
    objects.emplace_back(pool.Acquire());
  }

  ASSERT_EQ(10, pool.GetNumAllocatedObjects());
  ASSERT_EQ(0, pool.GetNumCachedObjects());
}

}  // namespace uom
