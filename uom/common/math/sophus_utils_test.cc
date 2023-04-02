// Copyright @2023 UOM project. All rights reserved.

#include "uom/common/math/sophus_utils.h"
#include <vector>

#include "gtest/gtest.h"

#include "uom/common/testing/eigen.h"
#include "uom/common/testing/jacobian.h"

TEST(SophusUtilsTest, Se3LogExp) {
  // All zeros
  {
    Sophus::SE3d se3;
    const auto actual = Se3Log(se3);
    EXPECT_EIGEN_EQ(Eigen::Vector3d::Zero(), actual.head<3>());
    EXPECT_EIGEN_EQ(Eigen::Vector3d::Zero(), actual.tail<3>());
  }

  // Only translation
  {
    const Eigen::Vector3d expect = Eigen::Vector3d(1, 2, 3);
    Sophus::SE3d se3(Sophus::SO3d::exp(Eigen::Vector3d::Zero()), expect);
    const auto actual = Se3Log(se3);
    EXPECT_EIGEN_EQ(expect, actual.head<3>());
    EXPECT_EIGEN_EQ(Eigen::Vector3d::Zero(), actual.tail<3>());
  }

  using Vector6d = Eigen::Matrix<double, 6, 1>;

  for (int i = 0; i < 100; ++i) {
    const Vector6d expect = Vector6d::Random();
    const Sophus::SE3d se3 = Se3Exp(expect);
    const Vector6d actual = Se3Log(se3);
    EXPECT_EIGEN_APPROX_EQ(expect, actual, 1e-6);
  }
}

TEST(SophusUtilsTest, LeftJacobianSO3) {
  const std::vector<Eigen::Vector3d> phi_list{
      Eigen::Vector3d(1e-13, 1e-13, 1e-13),
      Eigen::Vector3d(-1.2, 0.4, 0.9),
  };

  for (const auto& phi : phi_list) {
    const Sophus::SO3d inverse_exp_x = Sophus::SO3d::exp(phi).inverse();
    auto func = [&inverse_exp_x](const Eigen::Vector3d& phi_plus_delta) {
      return (Sophus::SO3d::exp(phi_plus_delta) * inverse_exp_x).log();
    };
    ASSERT_JACOBIAN_EQ(LeftJacobianSO3(phi), func, phi);
  }
}

TEST(SophusUtilsTest, LeftJacobianInverseSO3) {
  const std::vector<Eigen::Vector3d> phi_list{
      Eigen::Vector3d(1e-13, 1e-13, 1e-13),
      Eigen::Vector3d(-1.2, 0.4, 0.9),
  };

  for (const auto& phi : phi_list) {
    const Sophus::SO3d inverse_exp_x = Sophus::SO3d::exp(phi).inverse();
    auto func = [&inverse_exp_x](const Eigen::Vector3d& phi_plus_delta) {
      return (Sophus::SO3d::exp(phi_plus_delta) * inverse_exp_x).log();
    };
    ASSERT_JACOBIAN_EQ(LeftJacobianInverseSO3(phi).inverse(), func, phi);
  }
}

TEST(SophusUtilsTest, RightJacobianSO3) {
  const std::vector<Eigen::Vector3d> phi_list{
      Eigen::Vector3d(1e-13, 1e-13, 1e-13),
      Eigen::Vector3d(-1.2, 0.4, 0.9),
  };

  for (const auto& phi : phi_list) {
    const Sophus::SO3d inverse_exp_x = Sophus::SO3d::exp(phi).inverse();
    auto func = [&inverse_exp_x](const Eigen::Vector3d& phi_plus_delta) {
      return (inverse_exp_x * Sophus::SO3d::exp(phi_plus_delta)).log();
    };
    ASSERT_JACOBIAN_EQ(RightJacobianSO3(phi), func, phi);
  }
}

TEST(SophusUtilsTest, RightJacobianInverseSO3) {
  const std::vector<Eigen::Vector3d> phi_list{
      Eigen::Vector3d(1e-13, 1e-13, 1e-13),
      Eigen::Vector3d(-1.2, 0.4, 0.9),
  };

  for (const auto& phi : phi_list) {
    const Sophus::SO3d inverse_exp_x = Sophus::SO3d::exp(phi).inverse();
    auto func = [&inverse_exp_x](const Eigen::Vector3d& phi_plus_delta) {
      return (inverse_exp_x * Sophus::SO3d::exp(phi_plus_delta)).log();
    };
    ASSERT_JACOBIAN_EQ(RightJacobianInverseSO3(phi).inverse(), func, phi);
  }
}

