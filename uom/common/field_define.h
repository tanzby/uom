// Copyright @2023 UOM project. All rights reserved.

#pragma once

#define CONST_FIELD(type, name)                \
 public:                                       \
  const type& name() const { return name##_; } \
                                               \
 private:                                      \
  type name##_

#define FIELD(type, name)                                   \
 public:                                                    \
  type* mutable_##name() { return &name##_; }               \
  void set_##name(type name) { name##_ = std::move(name); } \
  CONST_FIELD(type, name)
