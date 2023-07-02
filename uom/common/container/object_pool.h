// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "uom/common/defines.h"

namespace uom {

template <typename T>
struct OnAquireDoNothing {
  void operator()(T*) const {}
};

template <typename T, typename OnAquire = OnAquireDoNothing<T>>
class ObjectPool final {
 public:
  using Ptr = std::unique_ptr<T, std::function<void(T*)>>;

  ObjectPool() = default;
  ~ObjectPool() = default;

  void Reserve(int capacity) {
    pool_.reserve(capacity);
    const int old_size = pool_.size() + num_allocated_objects_;
    const int require_size = capacity - old_size;
    if (require_size > 0) {
      for (int i = 0; i < require_size; ++i) {
        pool_.emplace_back(std::make_unique<T>());
      }
      num_allocated_objects_ += require_size;
    }
  }

  template <typename... Args>
  Ptr Acquire(Args&&... args) {
    Ptr object;
    if (pool_.empty()) {
      object = MakePtr(new T(std::forward<Args>(args)...));  // NOLINT
      ++num_allocated_objects_;
      return object;
    }
    object = MakePtr(pool_.back().release());
    pool_.pop_back();
    OnAquire()(object.get());
    return object;
  }

  [[nodiscard]] int GetNumCachedObjects() const { return pool_.size(); }
  [[nodiscard]] int GetNumAllocatedObjects() const { return num_allocated_objects_; }

 private:
  Ptr MakePtr(T* resource) {
    return Ptr(resource,
               std::bind(&ObjectPool::RecycleObjectCallback, this, std::placeholders::_1));
  }

  void RecycleObjectCallback(T* object_to_delete) {
    if (object_to_delete != nullptr) {
      pool_.emplace_back(object_to_delete);
    }
  }

  int num_allocated_objects_ = 0;

  // Memory owned. Use the default deleter.
  std::vector<std::unique_ptr<T>> pool_;

  DISALLOW_COPY_MOVE_AND_ASSIGN(ObjectPool);
};

}  // namespace uom
