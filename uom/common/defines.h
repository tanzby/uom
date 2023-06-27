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

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&) = delete;      \
  void operator=(const TypeName&) = delete

#define DISALLOW_COPY_MOVE_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&) = delete;           \
  void operator=(const TypeName&) = delete;     \
  TypeName(TypeName&&) = delete;                \
  void operator=(TypeName&&) = delete

// Compiler attributes
#if (defined(__GNUC__) || defined(__APPLE__)) && !defined(SWIG)
// Compiler supports GCC-style attributes
#define MUST_USE_RESULT __attribute__((warn_unused_result))
#else
// Non-GCC equivalents
#define MUST_USE_RESULT
#endif
