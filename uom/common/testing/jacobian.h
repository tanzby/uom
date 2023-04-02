// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <utility>
#include "Eigen/Eigen"
#include "gtest/gtest.h"

#define EXPECT_JACOBIAN_EQ(expect, func, x0) \
  EXPECT_PRED_FORMAT3(::testing::DefaultEqualsJacobian, expect, func, x0)
#define ASSERT_JACOBIAN_EQ(expect, func, x0) \
  ASSERT_PRED_FORMAT3(::testing::DefaultEqualsJacobian, expect, func, x0)

#define EXPECT_JACOBIAN_APPROX_EQ(expect, func, x0, eps, max_norm) \
  EXPECT_PRED_FORMAT5(::testing::ApproxEqualsJacobian, expect, func, x0, eps, max_norm)
#define ASSERT_JACOBIAN_APPROX_EQ(expect, func, x0, eps, max_norm) \
  ASSERT_PRED_FORMAT5(::testing::ApproxEqualsJacobian, expect, func, x0, eps, max_norm)

namespace testing {

namespace internal {
template <class Scalar>
struct TestConstants;

template <>
struct TestConstants<double> {
  static constexpr double epsilon = 1e-8;
  static constexpr double max_norm = 1e-3;
};

template <>
struct TestConstants<float> {
  static constexpr double epsilon = 1e-2;
  static constexpr double max_norm = 1e-2;
};
}  // namespace internal

template <typename Derived1, typename Derived2, typename F>
AssertionResult ApproxEqualsJacobian(const char* /*expect*/,
                                     const char* /*func*/,
                                     const char* /*x0*/,
                                     const char* /*eps*/,
                                     const char* /*max_norm*/,
                                     const Eigen::MatrixBase<Derived1>& expect,
                                     F&& func,
                                     const Eigen::MatrixBase<Derived2>& x0,
                                     double eps,
                                     double max_norm) {
  using Scalar = typename Derived1::Scalar;

  // Initialize the size of numerical_jacobian and set to zero.
  Eigen::MatrixX<Scalar> numerical_jacobian = expect;
  numerical_jacobian.setZero();

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> increment = x0;

  for (int i = 0; i < numerical_jacobian.cols(); ++i) {
    increment.setZero();
    increment[i] += eps;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> fpe = func(x0 + increment);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> fme = func(x0 - increment);
    numerical_jacobian.col(i) = (fpe - fme) / (2 * eps);
  }

  const Scalar diff = (expect - numerical_jacobian).template lpNorm<Eigen::Infinity>();
  if (diff <= max_norm) {
    return AssertionSuccess();
  }
  return AssertionFailure() << "Numerical jacobian is different from expect, max norm = " << diff
                            << " > " << max_norm << std::endl
                            << "expect:\n"
                            << expect << std::endl
                            << "numerical:\n"
                            << numerical_jacobian;
}

template <typename Derived1, typename Derived2, typename F>
AssertionResult DefaultEqualsJacobian(const char* /*expect*/,
                                      const char* /*func*/,
                                      const char* /*x0*/,
                                      const Eigen::MatrixBase<Derived1>& expect,
                                      F&& func,
                                      const Eigen::MatrixBase<Derived2>& x0) {
  return ApproxEqualsJacobian("",
                              "",
                              "",
                              "",
                              "",
                              expect,
                              std::forward<decltype(func)>(func),
                              x0,
                              internal::TestConstants<typename Derived1::Scalar>::epsilon,
                              internal::TestConstants<typename Derived1::Scalar>::max_norm);
}

}  // namespace testing
