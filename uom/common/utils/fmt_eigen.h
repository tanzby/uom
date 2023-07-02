// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include "Eigen/Eigen"
#include "fmt/ostream.h"

template <typename T>
struct fmt::formatter<T, std::enable_if_t<std::is_base_of_v<Eigen::DenseBase<T>, T>, char>>
    : ostream_formatter {};
