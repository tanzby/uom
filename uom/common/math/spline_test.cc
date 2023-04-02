// Copyright @2023 UOM project. All rights reserved.

#include "uom/common/math/spline.h"

#include <cstdint>

#include "gtest/gtest.h"

#include "uom/common/math/sophus_utils.h"
#include "uom/common/testing/eigen.h"
#include "uom/common/testing/jacobian.h"

namespace {

template <int DIM, int N, int DERIV>
void TestRdSplineEvaluate(const RdSpline<DIM, N>& spline, int64_t t_ns) {
  using VectorD = typename RdSpline<DIM, N>::VectorD;
  using MatrixD = typename RdSpline<DIM, N>::MatrixD;

  typename RdSpline<DIM, N>::JacobianStruct J_spline;

  spline.template Evaluate<DERIV>(t_ns, &J_spline);

  VectorD x0;
  x0.setZero();

  for (int64_t i = 0; i < 3 * N; i++) {
    std::stringstream ss;

    ss << "d_val_d_knot" << i << " time " << t_ns;

    MatrixD J_a;
    J_a.setZero();

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a.diagonal().setConstant(J_spline.d_val_d_knot[i - J_spline.start_idx]);
    }

    auto func = [&](const VectorD& x) {
      RdSpline<DIM, N> spline1 = spline;
      spline1.GetKnot(i) += x;
      return spline1.template Evaluate<DERIV>(t_ns);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0);
  }
}

template <int DIM, int N, int DERIV>
void TestRdSplineTimeDerivative(const RdSpline<DIM, N>& spline, int64_t t_ns) {
  using VectorD = typename RdSpline<DIM, N>::VectorD;

  const VectorD d_val_d_t = spline.template Evaluate<DERIV + 1>(t_ns);

  Eigen::Matrix<double, 1, 1> x0;
  x0.setZero();

  auto func = [&](const Eigen::Matrix<double, 1, 1>& x) {
    int64_t inc = x[0] * 1e9;
    return spline.template Evaluate<DERIV>(t_ns + inc);
  };
  EXPECT_JACOBIAN_EQ(d_val_d_t, func, x0);
}

template <int N>
void TestSo3SplineEvaluate(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<N>::Vector3;
  using MatrixD = typename So3Spline<N>::Matrix3;
  using SO3 = typename So3Spline<N>::SO3;

  typename So3Spline<N>::JacobianStruct J_spline;

  SO3 res = spline.Evaluate(t_ns, &J_spline);

  VectorD x0;
  x0.setZero();

  for (int i = 0; i < 3 * N; i++) {
    std::stringstream ss;

    MatrixD J_a;
    J_a.setZero();

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    }

    auto func = [&](const VectorD& x) {
      So3Spline<N> spline1 = spline;
      spline1.GetKnot(i) = SO3::exp(x) * spline.GetKnot(i);
      const SO3 res1 = spline1.Evaluate(t_ns);
      return (res1 * res.inverse()).log();
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "d_val_d_knot" << i << " time " << t_ns;
  }
}

template <int N>
void TestSo3SplineVelocity(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<N>::Vector3;
  using SO3 = typename So3Spline<N>::SO3;

  SO3 res = spline.Evaluate(t_ns);

  VectorD d_res_d_t = spline.VelocityBody(t_ns);

  Eigen::Matrix<double, 1, 1> x0;
  x0.setZero();

  auto func = [&](const Eigen::Matrix<double, 1, 1>& x) {
    int64_t inc = x[0] * 1e9;
    return (res.inverse() * spline.Evaluate(t_ns + inc)).log();
  };
  EXPECT_JACOBIAN_EQ(d_res_d_t, func, x0);
}

template <int N>
void TestSo3SplineAccelation(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<5>::Vector3;

  VectorD vel1;
  VectorD d_res_d_t = spline.AccelerationBody(t_ns, &vel1);

  VectorD vel2 = spline.VelocityBody(t_ns);
  EXPECT_TRUE(vel1.isApprox(vel2));

  Eigen::Matrix<double, 1, 1> x0;
  x0.setZero();

  auto func = [&](const Eigen::Matrix<double, 1, 1>& x) {
    int64_t inc = x[0] * 1e9;
    return spline.VelocityBody(t_ns + inc);
  };
  EXPECT_JACOBIAN_EQ(d_res_d_t, func, x0);
}

template <int N>
void TestSo3SplineJerk(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<5>::Vector3;

  VectorD vel1;
  VectorD accel1;
  VectorD d_res_d_t = spline.JerkBody(t_ns, &vel1, &accel1);

  VectorD vel2 = spline.VelocityBody(t_ns);
  EXPECT_TRUE(vel1.isApprox(vel2));

  VectorD accel2 = spline.AccelerationBody(t_ns);
  EXPECT_TRUE(accel1.isApprox(accel2));

  Eigen::Matrix<double, 1, 1> x0;
  x0.setZero();

  auto func = [&](const Eigen::Matrix<double, 1, 1>& x) {
    int64_t inc = x[0] * 1e9;
    return spline.AccelerationBody(t_ns + inc);
  };
  EXPECT_JACOBIAN_EQ(d_res_d_t, func, x0);
}

template <int N>
void TestSo3SplineEvaluateVel(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<5>::Vector3;
  using MatrixD = typename So3Spline<5>::Matrix3;
  using SO3 = typename So3Spline<5>::SO3;

  typename So3Spline<N>::JacobianStruct J_spline;

  VectorD res = spline.VelocityBody(t_ns, &J_spline);
  VectorD res_ref = spline.VelocityBody(t_ns);

  ASSERT_TRUE(res_ref.isApprox(res)) << "res and res_ref are not the same";

  VectorD x0;
  x0.setZero();

  for (int i = 0; i < 3 * N; i++) {
    std::stringstream ss;

    MatrixD J_a;
    J_a.setZero();

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    }

    auto func = [&](const VectorD& x) {
      So3Spline<N> spline1 = spline;
      spline1.GetKnot(i) = SO3::exp(x) * spline.GetKnot(i);

      return spline1.VelocityBody(t_ns);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "d_vel_d_knot" << i << " time " << t_ns;
  }
}

template <int N>
void TestSo3SplineEvaluateAccel(const So3Spline<N>& spline, int64_t t_ns) {
  using VectorD = typename So3Spline<N>::Vector3;
  using MatrixD = typename So3Spline<N>::Matrix3;
  using SO3 = typename So3Spline<N>::SO3;

  typename So3Spline<N>::JacobianStruct J_accel;
  typename So3Spline<N>::JacobianStruct J_vel;

  VectorD vel;
  VectorD vel_ref;

  VectorD res = spline.AccelerationBody(t_ns, &J_accel, &vel, &J_vel);
  VectorD res_ref = spline.AccelerationBody(t_ns, &vel_ref);

  ASSERT_TRUE(vel_ref.isApprox(vel)) << "vel and vel_ref are not the same";
  ASSERT_TRUE(res_ref.isApprox(res)) << "res and res_ref are not the same";

  VectorD x0;
  x0.setZero();

  // Test velocity Jacobian
  for (int i = 0; i < 3 * N; i++) {
    std::stringstream ss;

    MatrixD J_a;
    J_a.setZero();

    if (i >= J_vel.start_idx && i < J_vel.start_idx + N) {
      J_a = J_vel.d_val_d_knot[i - J_vel.start_idx];
    }

    auto func = [&](const VectorD& x) {
      So3Spline<N> spline1 = spline;
      spline1.GetKnot(i) = SO3::exp(x) * spline.GetKnot(i);

      return spline1.VelocityBody(t_ns);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "d_vel_d_knot" << i << " time " << t_ns;
  }

  // Test acceleration Jacobian
  for (int i = 0; i < 3 * N; i++) {
    std::stringstream ss;

    MatrixD J_a;
    J_a.setZero();

    if (i >= J_accel.start_idx && i < J_accel.start_idx + N) {
      J_a = J_accel.d_val_d_knot[i - J_accel.start_idx];
    }

    auto func = [&](const VectorD& x) {
      So3Spline<N> spline1 = spline;
      spline1.GetKnot(i) = SO3::exp(x) * spline.GetKnot(i);

      return spline1.AccelerationBody(t_ns);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "d_accel_d_knot" << i << " time " << t_ns;
  }
}

template <int N>
void TestSe3SplineGyroResidual(const Se3Spline<N>& s, int64_t t_ns) {
  typename Se3Spline<N>::SO3JacobianStruct J_spline;
  Eigen::Matrix<double, 3, 12> J_bias;

  CalibGyroBias<double> bias;

  bias.SetRandom();

  // bias << 0.01, -0.02, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  Eigen::Vector3d measurement = s.RotationalVelocityBody(t_ns);

  s.GyroResidual(t_ns, measurement, bias, &J_spline, &J_bias);

  for (int64_t i = 0; i < s.NumKnots(); i++) {
    Sophus::Vector3d x0;
    x0.setZero();

    Sophus::Matrix3d J_a;

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    } else {
      J_a.setZero();
    }

    auto func = [&](const Sophus::Vector3d& x) {
      Se3Spline<N> s1 = s;
      s1.GetKnotSO3(i) = Sophus::SO3d::exp(x) * s.GetKnotSO3(i);
      return s1.GyroResidual(t_ns, measurement, bias);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0)
        << "Spline order " << N << " d_gyro_res_d_knot" << i << " time " << t_ns;
  }

  {
    Eigen::Matrix<double, 12, 1> x0;
    x0.setZero();

    auto func = [&](const Eigen::Matrix<double, 12, 1>& x) {
      auto b1 = bias;
      b1 += x;
      return s.GyroResidual(t_ns, measurement, b1);
    };
    EXPECT_JACOBIAN_EQ(J_bias, func, x0) << "Spline order " << N << " d_gyro_res_d_bias";
  }
}

template <int N>
void TestSe3SplineAccelResidual(const Se3Spline<N>& s, int64_t t_ns) {
  typename Se3Spline<N>::AccelPosSO3JacobianStruct J_spline;
  Eigen::Matrix3d J_g;

  Eigen::Matrix<double, 3, 9> J_bias;

  CalibAccelBias<double> bias;
  bias.SetRandom();
  // bias << -0.4, 0.5, -0.6, 0, 0, 0, 0, 0, 0;

  Eigen::Vector3d g(0, 0, -9.81);
  Eigen::Vector3d measurement = s.LinearAccelerationWorld(t_ns) + g;

  s.AccelResidual(t_ns, measurement, bias, g, &J_spline, &J_bias, &J_g);

  for (int64_t i = 0; i < s.NumKnots(); i++) {
    Sophus::Vector6d x0;
    x0.setZero();

    typename Se3Spline<N>::Mat36 J_a;

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    } else {
      J_a.setZero();
    }

    auto func = [&](const Sophus::Vector6d& x) {
      Se3Spline<N> s1 = s;
      s1.BoxPlusForKnot(i, x);

      return s1.AccelResidual(t_ns, measurement, bias, g);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "Spline order " << N << " d_accel_res_d_knot" << i;
  }

  {
    Eigen::Matrix<double, 9, 1> x0;
    x0.setZero();

    auto func = [&](const Eigen::Matrix<double, 9, 1>& x) {
      auto b1 = bias;
      b1 += x;
      return s.AccelResidual(t_ns, measurement, b1, g);
    };
    EXPECT_JACOBIAN_EQ(J_bias, func, x0) << "Spline order " << N << " d_accel_res_d_bias";
  }

  {
    Sophus::Vector3d x0;
    x0.setZero();

    auto func = [&](const Sophus::Vector3d& x) {
      return s.AccelResidual(t_ns, measurement, bias, g + x);
    };
    EXPECT_JACOBIAN_EQ(J_g, func, x0) << "Spline order " << N << " d_accel_res_d_g";
  }
}

template <int N>
void TestSe3SplineOrientationResidual(const Se3Spline<N>& s, int64_t t_ns) {
  typename Se3Spline<N>::SO3JacobianStruct J_spline;

  Sophus::SO3d measurement = s.Pose(t_ns).so3();

  s.OrientationResidual(t_ns, measurement, &J_spline);

  for (int64_t i = 0; i < s.NumKnots(); i++) {
    Sophus::Matrix3d J_a;

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    } else {
      J_a.setZero();
    }

    Sophus::Vector3d x0;
    x0.setZero();

    auto func = [&](const Sophus::Vector3d& x_rot) {
      Sophus::Vector6d x;
      x.setZero();
      x.tail<3>() = x_rot;
      Se3Spline<N> s1 = s;
      s1.BoxPlusForKnot(i, x);
      return s1.OrientationResidual(t_ns, measurement);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "Spline order " << N << " d_rot_res_d_knot" << i;
  }
}

template <int N>
void TestSe3SplinePositionResidual(const Se3Spline<N>& s, int64_t t_ns) {
  typename Se3Spline<N>::PosJacobianStruct J_spline;

  Eigen::Vector3d measurement = s.Pose(t_ns).translation();

  s.PositionResidual(t_ns, measurement, &J_spline);

  for (int64_t i = 0; i < s.NumKnots(); i++) {
    Sophus::Matrix3d J_a;
    J_a.setZero();

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a.diagonal().setConstant(J_spline.d_val_d_knot[i - J_spline.start_idx]);
    }

    Sophus::Vector3d x0;
    x0.setZero();

    auto func = [&](const Sophus::Vector3d& x_rot) {
      Sophus::Vector6d x;
      x.setZero();
      x.head<3>() = x_rot;

      Se3Spline<N> s1 = s;
      s1.BoxPlusForKnot(i, x);

      return s1.PositionResidual(t_ns, measurement);
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "Spline order " << N << " d_pos_res_d_knot" << i;
  }
}

template <int N>
void TestSe3SplinePose(const Se3Spline<N>& s, int64_t t_ns) {
  typename Se3Spline<N>::PosePosSO3JacobianStruct J_spline;

  Sophus::SE3d res = s.Pose(t_ns, &J_spline);

  Sophus::Vector6d x0;
  x0.setZero();

  for (int64_t i = 0; i < s.NumKnots(); i++) {
    typename Se3Spline<N>::Mat6 J_a;

    if (i >= J_spline.start_idx && i < J_spline.start_idx + N) {
      J_a = J_spline.d_val_d_knot[i - J_spline.start_idx];
    } else {
      J_a.setZero();
    }

    auto func = [&](const Sophus::Vector6d& x) {
      Se3Spline<N> s1 = s;
      s1.BoxPlusForKnot(i, x);
      return Se3Log(res.inverse() * s1.Pose(t_ns));
    };
    EXPECT_JACOBIAN_EQ(J_a, func, x0) << "Spline order " << N << " d_pose_d_knot" << i;
  }

  {
    Eigen::Matrix<double, 1, 1> x0;
    x0[0] = 0;

    typename Se3Spline<N>::Vec6 J_pose_time;

    s.d_pose_d_t(t_ns, J_pose_time);

    auto func = [&](const Eigen::Matrix<double, 1, 1>& x) {
      int64_t t_ns_new = t_ns;
      t_ns_new += x[0] * 1e9;
      return Se3Log(res.inverse() * s.Pose(t_ns_new));
    };
    EXPECT_JACOBIAN_EQ(J_pose_time, func, x0);
  }
}

}  // namespace

TEST(So3SplineTest, ComputeBaseCoefficients) {
  Eigen::Matrix<double, 5, 5> expect;
  // clang-format off
  expect << 
  1, 1, 1, 1, 1,
  0, 1, 2, 3, 4,
  0, 0, 2, 6, 12,
  0, 0, 0, 6, 24,
  0, 0, 0, 0, 24;
  // clang-format on

  const auto actual = ComputeBaseCoefficients<5>();
  EXPECT_EIGEN_EQ(expect, actual);
}

TEST(RdSplineTest, EvaluateKnots4) {
  static constexpr int DIM = 3;
  static constexpr int N = 4;
  static constexpr int Deriv = 0;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnots5) {
  static constexpr int DIM = 3;
  static constexpr int N = 5;
  static constexpr int Deriv = 0;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnots6) {
  static constexpr int DIM = 3;
  static constexpr int N = 6;
  static constexpr int Deriv = 0;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsVelocity4) {
  static constexpr int DIM = 3;
  static constexpr int N = 4;
  static constexpr int Deriv = 1;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsVelocity5) {
  static constexpr int DIM = 3;
  static constexpr int N = 5;
  static constexpr int Deriv = 1;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsVelocity6) {
  static constexpr int DIM = 3;
  static constexpr int N = 6;
  static constexpr int Deriv = 1;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsAcceleration4) {
  static constexpr int DIM = 3;
  static constexpr int N = 4;
  static constexpr int Deriv = 2;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsAcceleration5) {
  static constexpr int DIM = 3;
  static constexpr int N = 5;
  static constexpr int Deriv = 2;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateKnotsAcceleration6) {
  static constexpr int DIM = 3;
  static constexpr int N = 6;
  static constexpr int Deriv = 2;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestRdSplineEvaluate<DIM, N, Deriv>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateTimeDerivative4) {
  static constexpr int DIM = 3;
  static constexpr int N = 4;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 0>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateTimeDerivative5) {
  static constexpr int DIM = 3;
  static constexpr int N = 5;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 0>(spline, t_ns);
  }
}

TEST(RdSplineTest, EvaluateTimeDerivative6) {
  static constexpr int DIM = 3;
  static constexpr int N = 6;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 0>(spline, t_ns);
  }
}

TEST(RdSplineTest, VelocityTimeDerivative4) {
  static constexpr int DIM = 3;
  static constexpr int N = 4;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 1>(spline, t_ns);
  }
}

TEST(RdSplineTest, VelocityTimeDerivative5) {
  static constexpr int DIM = 3;
  static constexpr int N = 5;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 1>(spline, t_ns);
  }
}

TEST(RdSplineTest, VelocityTimeDerivative6) {
  static constexpr int DIM = 3;
  static constexpr int N = 6;

  RdSpline<DIM, N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestRdSplineTimeDerivative<DIM, N, 1>(spline, t_ns);
  }
}

TEST(So3SplineTest, EvaluateKnots4) {
  static constexpr int N = 4;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluate(spline, t_ns);
  }
}

TEST(So3SplineTest, EvaluateKnots5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluate(spline, t_ns);
  }
}

TEST(So3SplineTest, EvaluateKnots6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluate(spline, t_ns);
  }
}

TEST(So3SplineTest, Velocity4) {
  static constexpr int N = 4;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineVelocity(spline, t_ns);
  }
}

TEST(So3SplineTest, Velocity5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineVelocity(spline, t_ns);
  }
}

TEST(So3SplineTest, Velocity6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineVelocity(spline, t_ns);
  }
}

TEST(So3SplineTest, Acceleration4) {
  static constexpr int N = 4;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineAccelation(spline, t_ns);
  }
}

TEST(So3SplineTest, Acceleration5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineAccelation(spline, t_ns);
  }
}

TEST(So3SplineTest, Acceleration6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineAccelation(spline, t_ns);
  }
}

TEST(So3SplineTest, Jerk5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineJerk(spline, t_ns);
  }
}

TEST(So3SplineTest, Jerk6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 1e8; t_ns < spline.MaxTimeNs() - 1e8; t_ns += 1e8) {
    TestSo3SplineJerk(spline, t_ns);
  }
}

TEST(So3SplineTest, VelocityKnots4) {
  static constexpr int N = 4;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateVel(spline, t_ns);
  }
}

TEST(So3SplineTest, VelocityKnots5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateVel(spline, t_ns);
  }
}

TEST(So3SplineTest, VelocityKnots6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateVel(spline, t_ns);
  }
}

TEST(So3SplineTest, AccelerationKnots4) {
  static constexpr int N = 4;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateAccel(spline, t_ns);
  }
}

TEST(So3SplineTest, AccelerationKnots5) {
  static constexpr int N = 5;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateAccel(spline, t_ns);
  }
}

TEST(So3SplineTest, AccelerationKnots6) {
  static constexpr int N = 6;

  So3Spline<N> spline(int64_t(2e9));
  spline.GenerateRandomTrajectory(3 * N);

  for (int64_t t_ns = 0; t_ns < spline.MaxTimeNs(); t_ns += 1e8) {
    TestSo3SplineEvaluateAccel(spline, t_ns);
  }
}

TEST(Se3SplineTest, GyroResidual) {
  static constexpr int N = 5;

  const int num_knots = 3 * N;
  Se3Spline<N> s(int64_t(2e9));
  s.GenerateRandomTrajectory(num_knots);

  for (int64_t t_ns = 0; t_ns < s.MaxTimeNs(); t_ns += 1e8) {
    TestSe3SplineGyroResidual<N>(s, t_ns);
  }
}

TEST(Se3SplineTest, AccelResidual) {
  static constexpr int N = 5;

  const int num_knots = 3 * N;
  Se3Spline<N> s(int64_t(2e9));
  s.GenerateRandomTrajectory(num_knots);

  for (int64_t t_ns = 0; t_ns < s.MaxTimeNs(); t_ns += 1e8) {
    TestSe3SplineAccelResidual<N>(s, t_ns);
  }
}

TEST(Se3SplineTest, PositionResidual) {
  static constexpr int N = 5;

  const int num_knots = 3 * N;
  Se3Spline<N> s(int64_t(2e9));
  s.GenerateRandomTrajectory(num_knots);

  for (int64_t t_ns = 0; t_ns < s.MaxTimeNs(); t_ns += 1e8) {
    TestSe3SplinePositionResidual<N>(s, t_ns);
  }
}

TEST(Se3SplineTest, OrientationResidual) {
  static constexpr int N = 5;

  const int num_knots = 3 * N;
  Se3Spline<N> s(int64_t(2e9));
  s.GenerateRandomTrajectory(num_knots);

  for (int64_t t_ns = 0; t_ns < s.MaxTimeNs(); t_ns += 1e8) {
    TestSe3SplineOrientationResidual<N>(s, t_ns);
  }
}

TEST(Se3SplineTest, Pose) {
  static constexpr int N = 5;

  const int num_knots = 3 * N;
  Se3Spline<N> s(int64_t(2e9));
  s.GenerateRandomTrajectory(num_knots);

  int64_t offset = 100;

  for (int64_t t_ns = offset; t_ns < s.MaxTimeNs() - offset; t_ns += 1e8) {
    TestSe3SplinePose<N>(s, t_ns);
  }
}
