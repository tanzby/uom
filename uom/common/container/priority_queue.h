// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

namespace uom {

// Keeps font element is always the smallest element.
template <typename Movable, typename Less = std::less<Movable>>
class PriorityQueue final {
 public:
  // Returns size of queue.
  [[nodiscard]] int size() const noexcept {
    return static_cast<int>(elements_.end() - elements_.begin());
  }

  // Test whether container is empty.
  [[nodiscard]] bool empty() const noexcept { return elements_.empty(); }

  // Accesses top element.
  const Movable& top() const { return elements_.front(); }

  // Clears the queue.
  void clear() { elements_.clear(); }

  // Reserves memory space for container.
  void reserve(int64_t capcity) { elements_.reserve(capcity); }

  // Removes top element.
  Movable pop() {
    std::pop_heap(elements_.begin(), elements_.end(), greater_);
    Movable result = std::move(elements_.back());
    elements_.pop_back();
    return std::move(result);
  }

  // Inserts element.
  void push(Movable element) {
    elements_.emplace_back(std::move(element));
    std::push_heap(elements_.begin(), elements_.end(), greater_);
  }

  // Constructs and insert element.
  template <typename... Args>
  void emplace(Args&&... args) {
    elements_.emplace_back(std::forward<Args>(args)...);
    std::push_heap(elements_.begin(), elements_.end(), greater_);
  }

  // Swaps containers.
  void swap(PriorityQueue& rhs) noexcept { elements_.swap(rhs); }

 private:
  struct Greater {
    inline bool operator()(const Movable& lhs, const Movable& rhs) const { return less(rhs, lhs); }
    Less less;
  };

  Greater greater_;
  std::vector<Movable> elements_;
};

}  // namespace uom
