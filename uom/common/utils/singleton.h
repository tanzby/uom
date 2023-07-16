// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include "uom/common/defines.h"

namespace uom {

template <typename T>
class Singleton {
 public:
  static T* GetInstance() {
    static T instance;
    return &instance;
  }

 protected:
  Singleton() = default;
  virtual ~Singleton() = default;

 private:
  DISALLOW_COPY_MOVE_AND_ASSIGN(Singleton);
};

}  // namespace uom
