// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <limits>
#include <type_traits>

#include "Eigen/Core"
#include "gtest/gtest.h"

#define EXPECT_EIGEN_EQ(expect, actual) EXPECT_PRED_FORMAT2(::testing::EqualsEigen, expect, actual)
#define ASSERT_EIGEN_EQ(expect, actual) ASSERT_PRED_FORMAT2(::testing::EqualsEigen, expect, actual)
#define EXPECT_EIGEN_APPROX_EQ(expect, actual, eps) \
  EXPECT_PRED_FORMAT3(::testing::ApproxEqualsEigen, expect, actual, eps)
#define ASSERT_EIGEN_APPROX_EQ(expect, actual, eps) \
  ASSERT_PRED_FORMAT3(::testing::ApproxEqualsEigen, expect, actual, eps)

namespace testing {

template <typename Derived1, typename Derived2>
AssertionResult ApproxEqualsEigen(const char* /*expect*/,
                                  const char* /*actual*/,
                                  const char* /*eps*/,
                                  const Eigen::MatrixBase<Derived1>& expect,
                                  const Eigen::MatrixBase<Derived2>& actual,
                                  double eps) {

  if (expect.rows() != actual.rows() || expect.cols() != actual.cols()) {
    return AssertionFailure()
           << "The expect and actual matrix isn't (approx) equal as they have different size: \n"
           << "expect: rows=" << expect.rows() << ", cols=" << expect.cols() << std::endl
           << "expect: rows=" << actual.rows() << ", cols=" << actual.cols() << std::endl;
  }

  const typename Derived1::Scalar diff = (expect - actual).template lpNorm<Eigen::Infinity>();

  if (diff <= eps) {
    return AssertionSuccess();
  }

  return AssertionFailure() << "The expect and actual matrix isn't (approx) equal, diff = " << diff
                            << " > " << eps << ", where"
                            << "\nexpect evaluates to:\n"
                            << expect << "\nactual evaluates to:\n"
                            << actual;
}

template <typename Derived1, typename Derived2>
AssertionResult EqualsEigen(const char* /*expect*/,
                            const char* /*actual*/,
                            const Eigen::MatrixBase<Derived1>& expect,
                            const Eigen::MatrixBase<Derived2>& actual) {
  return ApproxEqualsEigen(
      "", "", "", expect, actual, std::numeric_limits<typename Derived1::Scalar>::epsilon());
}

}  // namespace testing