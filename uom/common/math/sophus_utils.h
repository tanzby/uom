// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include "sophus/se3.hpp"
#include "glog/logging.h"

// Decoupled version of expmap for SE(3). Returns [tx, ty, tz, rx, ry, rz].
template <typename Scalar>
inline typename Sophus::SE3<Scalar>::Tangent Se3Log(const Sophus::SE3<Scalar>& s) {
  typename Sophus::SE3<Scalar>::Tangent xi;
  xi.template tail<3>() = s.so3().log();
  xi.template head<3>() = s.translation();
  return xi;
}

// Decoupled version of logmap for SE(3). Accepts [tx, ty, tz, rx, ry, rz].
template <typename Derived>
inline Sophus::SE3<typename Derived::Scalar> Se3Exp(const Eigen::MatrixBase<Derived>& xi) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 6);
  using Scalar = typename Derived::Scalar;
  return Sophus::SE3<Scalar>(Sophus::SO3<Scalar>::exp(xi.template tail<3>()),
                             xi.template head<3>());
}

// Returns the left Jacobian of SO(3).
//
//           log(Exp(phi+mu) * Exp(phi)^{-1})
// J_l = lim ---------------------------------
//       mu>0               mu
//
template <typename Derived>
inline Eigen::Matrix3<typename Derived::Scalar> LeftJacobianSO3(
    const Eigen::MatrixBase<Derived>& phi) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
  return Sophus::SO3<typename Derived::Scalar>::leftJacobian(phi);
}

// Returns the right Jacobian of SO(3).
//
//           log(Exp(phi)^{-1} * Exp(phi+mu))
// J_l = lim ---------------------------------
//       mu>0               mu
//
template <typename Derived>
inline Eigen::Matrix3<typename Derived::Scalar> RightJacobianSO3(
    const Eigen::MatrixBase<Derived>& phi) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);

  using Scalar = typename Derived::Scalar;

  const Scalar phi_norm2 = phi.squaredNorm();
  const Eigen::Matrix3<Scalar> phi_hat = Sophus::SO3<Scalar>::hat(phi);
  const Eigen::Matrix3<Scalar> phi_hat2 = phi_hat * phi_hat;

  Eigen::Matrix3<Scalar> J_phi = Eigen::Matrix3<Scalar>::Identity();

  if (phi_norm2 > Sophus::Constants<Scalar>::epsilon()) {
    Scalar phi_norm = std::sqrt(phi_norm2);
    Scalar phi_norm3 = phi_norm2 * phi_norm;
    J_phi -= phi_hat * (1 - std::cos(phi_norm)) / phi_norm2;
    J_phi += phi_hat2 * (phi_norm - std::sin(phi_norm)) / phi_norm3;
  } else {
    J_phi -= phi_hat / 2;
    J_phi += phi_hat2 / 6;
  }

  return J_phi;
}

// Returns the left Inverse Jacobian of SO(3).
template <typename Derived>
inline Eigen::Matrix3<typename Derived::Scalar> LeftJacobianInverseSO3(
    const Eigen::MatrixBase<Derived>& phi) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
  return Sophus::SO3<typename Derived::Scalar>::leftJacobianInverse(phi);
}

// Returns the right Inverse Jacobian of SO(3).
template <typename Derived>
inline Eigen::Matrix3<typename Derived::Scalar> RightJacobianInverseSO3(
    const Eigen::MatrixBase<Derived>& phi) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);

  using Scalar = typename Derived::Scalar;

  const Scalar phi_norm2 = phi.squaredNorm();
  const Eigen::Matrix3<Scalar> phi_hat = Sophus::SO3<Scalar>::hat(phi);
  const Eigen::Matrix3<Scalar> phi_hat2 = phi_hat * phi_hat;

  Eigen::Matrix3<Scalar> J_phi = Eigen::Matrix3<Scalar>::Identity();
  J_phi += phi_hat / 2;

  if (phi_norm2 > Sophus::Constants<Scalar>::epsilon()) {
    Scalar phi_norm = std::sqrt(phi_norm2);

    // We require that the angle is in range [0, pi]. We check if we are close
    // to pi and apply a Taylor expansion to scalar multiplier of phi_hat2.
    // Technically, log(exp(phi)exp(epsilon)) is not continuous / differentiable
    // at phi=pi, but we still aim to return a reasonable value for all valid
    // inputs.
    CHECK_LE(phi_norm, M_PI + Sophus::Constants<Scalar>::epsilon());

    if (phi_norm < M_PI - Sophus::Constants<Scalar>::epsilonSqrt()) {
      // regular case for range (0,pi)
      J_phi += phi_hat2 *
               (1 / phi_norm2 - (1 + std::cos(phi_norm)) / (2 * phi_norm * std::sin(phi_norm)));
    } else {
      // 0th-order Taylor expansion around pi
      J_phi += phi_hat2 / (M_PI * M_PI);
    }
  } else {
    // Taylor expansion around 0
    J_phi += phi_hat2 / 12;
  }

  return J_phi;
}

// rightJacobianInvSO3
// rightJacobianSO3