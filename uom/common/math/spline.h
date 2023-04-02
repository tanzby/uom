// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <deque>
#include <iostream>

#include "Eigen/Dense"
#include "fmt/format.h"
#include "glog/logging.h"
#include "sophus/so3.hpp"

#include "uom/common/field_define.h"
#include "uom/common/math/sophus_utils.h"
#include "uom/common/math/spline_common.h"
#include "uom/sensors/calib_bias.h"

template <int D, int N, typename Scalar = double>
class SplineBase {
 public:
  using MatrixN = Eigen::Matrix<Scalar, N, N>;
  using VectorN = Eigen::Matrix<Scalar, N, 1>;
  using VectorD = Eigen::Matrix<Scalar, D, 1>;
  using MatrixD = Eigen::Matrix<Scalar, D, D>;

  /// @brief Vector of derivatives of time polynomial.
  ///
  /// Computes a derivative of \f$ \begin{bmatrix}1 & t & t^2 & \dots &
  /// t^{N-1}\end{bmatrix} \f$ with repect to time. For example, the first
  /// derivative would be \f$ \begin{bmatrix}0 & 1 & 2 t & \dots & (N-1)
  /// t^{N-2}\end{bmatrix} \f$.
  /// @param Derivative derivative to evaluate
  /// @param[in] t
  /// @param[out] result vector to store the result
  template <int Derivative>
  static void BaseCoeffsWithTime(Scalar t, VectorN* result) {
    CHECK_NOTNULL(result)->setZero();

    if constexpr (Derivative < N) {
      (*result)[Derivative] = kBaseCoeffs(Derivative, Derivative);

      Scalar ti = t;
      for (int j = Derivative + 1; j < N; j++) {
        (*result)[j] = kBaseCoeffs(Derivative, j) * ti;
        ti = ti * t;
      }
    }
  }

  // Base coefficients matrix.
  inline static const MatrixN kBaseCoeffs = ComputeBaseCoefficients<N, Scalar>();

  // Blending matrix.
  inline static const MatrixN kBlendingMatrix = ComputeBlendingMatrix<N, Scalar, false>();

  // Blending matrix.
  inline static const MatrixN kCumulativeBlendingMatrix = ComputeBlendingMatrix<N, Scalar, true>();
};

/// @brief Uniform B-spline for euclidean vectors with dimention D of order N.
///
/// @param D Dimension of euclidean vector space.
/// @param N Order of the spline.
/// @param Scalar Underlying floating type of this spline.
template <int D, int N, typename Scalar = double>
class RdSpline : public SplineBase<D, N, Scalar> {
  using Base = SplineBase<D, N, Scalar>;

 public:
  static constexpr Scalar NS_TO_S = 1e-9;  // Nanosecond to second conversion.
  static constexpr Scalar S_TO_NS = 1e9;   // Second to nanosecond conversion.

  using MatrixN = typename Base::MatrixN;
  using VectorN = typename Base::VectorN;
  using VectorD = typename Base::VectorD;
  using MatrixD = typename Base::MatrixD;

  // Struct to store the Jacobian of the spline. Since B-spline of order N has local support (only N
  // knots infuence the value) the Jacobian is zero for all knots except maximum N for value and all
  // derivatives.
  struct JacobianStruct {
    int64_t start_idx;                   // Start index of the non-zero elements of the Jacobian.
    std::array<Scalar, N> d_val_d_knot;  // Value of non-zero Jacobians.
  };

  /// @brief Default constructor
  RdSpline() = default;

  /// @brief Constructor with knot interval and start time.
  ///
  /// @param[in] time_interval_ns knot time interval in nanoseconds
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  RdSpline(int64_t time_interval_ns, int64_t start_time_ns = 0)
      : start_t_ns_(start_time_ns), dt_ns_(time_interval_ns) {
    pow_inv_dt_[0] = 1.0;
    pow_inv_dt_[1] = S_TO_NS / dt_ns_;

    for (int64_t i = 2; i < N; i++) {
      pow_inv_dt_[i] = pow_inv_dt_[i - 1] * pow_inv_dt_[1];
    }
  }

  /// @brief Cast to different scalar type
  template <typename Scalar2>
  RdSpline<D, N, Scalar2> Cast() const {
    RdSpline<D, N, Scalar2> res;

    res.dt_ns_ = dt_ns_;
    res.start_t_ns_ = start_t_ns_;

    for (int i = 0; i < N; i++) {
      res.pow_inv_dt_[i] = pow_inv_dt_[i];
    }

    for (const auto& k : knots_) {
      res.knots_.emplace_back(k.template cast<Scalar2>());
    }

    return res;
  }

  /// @brief Maximum time represented by spline
  ///
  /// @return maximum time represented by spline in nanoseconds
  int64_t MaxTimeNs() const { return start_t_ns_ + (knots_.size() - N + 1) * dt_ns_ - 1; }

  /// @brief Minimum time represented by spline
  ///
  /// @return minimum time represented by spline in nanoseconds
  int64_t MinTimeNs() const { return start_t_ns_; }

  /// @brief Gererate random trajectory
  ///
  /// @param[in] n number of knots to generate
  /// @param[in] static_init if true the first N knots will be the same resulting in static initial
  /// condition
  void GenerateRandomTrajectory(int n, bool static_init = false) {
    if (static_init) {
      VectorD rnd = VectorD::Random() * 5;

      for (int i = 0; i < N; i++) {
        knots_.push_back(rnd);
      }
      for (int i = 0; i < n - N; i++) {
        knots_.push_back(VectorD::Random() * 5);
      }
    } else {
      for (int i = 0; i < n; i++) {
        knots_.push_back(VectorD::Random() * 5);
      }
    }
  }

  /// @brief Add knot to the end of the spline
  ///
  /// @param[in] knot knot to add
  void KnotsPushBack(const VectorD& knot) { knots_.push_back(knot); }

  /// @brief Remove knot from the back of the spline
  void KnotsPopBack() { knots_.pop_back(); }

  /// @brief Return the first knot of the spline
  ///
  /// @return first knot of the spline
  const VectorD& KnotsFront() const { return knots_.front(); }

  /// @brief Remove first knot of the spline and increase the start time
  void KnotsPopFront() {
    start_t_ns_ += dt_ns_;
    knots_.pop_front();
  }

  /// @brief Resize containter with knots
  ///
  /// @param[in] n number of knots
  void Resize(int64_t n) { knots_.resize(n); }

  /// @brief Return reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return reference to the knot
  VectorD& GetKnot(int i) { return knots_[i]; }

  /// @brief Return const reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return const reference to the knot
  const VectorD& GetKnot(int i) const { return knots_[i]; }

  /// @brief Evaluate value or derivative of the spline
  ///
  /// @param Derivative derivative to evaluate (0 for value)
  /// @param[in] time_ns time for evaluating of the spline in nanoseconds
  /// @param[out] J if not nullptr, return the Jacobian of the value with respect to knots
  /// @return value of the spline or derivative. Euclidean vector of dimention DIM.
  template <int Derivative = 0>
  VectorD Evaluate(int64_t time_ns, JacobianStruct* J = nullptr) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<Derivative>(u, &p);

    VectorN coeff = pow_inv_dt_[Derivative] * (this->kBlendingMatrix * p);

    VectorD res;
    res.setZero();

    for (int i = 0; i < N; i++) {
      res += coeff[i] * knots_[s + i];
      if (J != nullptr) {
        J->d_val_d_knot[i] = coeff[i];
      }
    }

    if (J != nullptr) {
      J->start_idx = s;
    }

    return res;
  }

  /// @brief Alias for first derivative of spline. See \ref evaluate.
  VectorD Velocity(int64_t time_ns, JacobianStruct* J = nullptr) const {
    return Evaluate<1>(time_ns, J);
  }

  /// @brief Alias for second derivative of spline. See \ref evaluate.
  VectorD Acceleration(int64_t time_ns, JacobianStruct* J = nullptr) const {
    return Evaluate<2>(time_ns, J);
  }

 private:
  using KnotQueue = std::deque<VectorD>;
  FIELD(int64_t, start_t_ns) = 0;     // Start time in nanoseconds.
  CONST_FIELD(int64_t, dt_ns) = 0;    // Knot interval in nanoseconds.
  CONST_FIELD(KnotQueue, knots);      // Knots.
  std::array<Scalar, N> pow_inv_dt_;  // Array with inverse powers of dt.
};

/// @brief Uniform cummulative B-spline for SO(3) of order N
///
/// @param N Order of the spline.
/// @param Scalar Underlying floating type of this spline.
template <int N, typename Scalar = double>
class So3Spline : public SplineBase<3, N, Scalar> {
  using Base = SplineBase<3, N, Scalar>;

 public:
  static constexpr int DEG = N - 1;  // Degree of the spline.

  static constexpr Scalar NS_TO_S = 1e-9;  // Nanosecond to second conversion
  static constexpr Scalar S_TO_NS = 1e9;   // Second to nanosecond conversion

  using MatrixN = typename Base::MatrixN;
  using VectorN = typename Base::VectorN;
  using Vector3 = typename Base::VectorD;
  using Matrix3 = typename Base::MatrixD;

  using SO3 = Sophus::SO3<Scalar>;

  /// @brief Struct to store the Jacobian of the spline
  ///
  /// Since B-spline of order N has local support (only N knots infuence the
  /// value) the Jacobian is zero for all knots except maximum N for value and
  /// all derivatives.
  struct JacobianStruct {
    int start_idx;
    std::array<Matrix3, N> d_val_d_knot;
  };

  /// @brief Constructor with knot interval and start time
  ///
  /// @param[in] time_interval_ns knot time interval in nanoseconds
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  So3Spline(int64_t time_interval_ns, int64_t start_time_ns = 0)
      : dt_ns_(time_interval_ns), start_t_ns_(start_time_ns) {
    pow_inv_dt_[0] = 1.0;
    pow_inv_dt_[1] = S_TO_NS / dt_ns_;
    pow_inv_dt_[2] = pow_inv_dt_[1] * pow_inv_dt_[1];
    pow_inv_dt_[3] = pow_inv_dt_[2] * pow_inv_dt_[1];
  }

  /// @brief Maximum time represented by spline
  ///
  /// @return maximum time represented by spline in nanoseconds
  int64_t MaxTimeNs() const { return start_t_ns_ + (knots_.size() - N + 1) * dt_ns_ - 1; }

  /// @brief Minimum time represented by spline
  ///
  /// @return minimum time represented by spline in nanoseconds
  int64_t MinTimeNs() const { return start_t_ns_; }

  /// @brief Gererate random trajectory
  ///
  /// @param[in] n number of knots to generate
  /// @param[in] static_init if true the first N knots will be the same
  /// resulting in static initial condition
  void GenerateRandomTrajectory(int n, bool static_init = false) {
    if (static_init) {
      Vector3 rnd = Vector3::Random() * M_PI;

      for (int i = 0; i < N; i++) {
        knots_.push_back(SO3::exp(rnd));
      }

      for (int i = 0; i < n - N; i++) {
        knots_.push_back(knots_.back() * SO3::exp(Vector3::Random() * M_PI / 2));
      }

    } else {
      knots_.push_back(SO3::exp(Vector3::Random() * M_PI));

      for (int i = 1; i < n; i++) {
        knots_.push_back(knots_.back() * SO3::exp(Vector3::Random() * M_PI / 2));
      }
    }
  }

  /// @brief Add knot to the end of the spline
  ///
  /// @param[in] knot knot to add
  void KnotsPushBack(const SO3& knot) { knots_.push_back(knot); }

  /// @brief Remove knot from the back of the spline
  void KnotsPopBack() { knots_.pop_back(); }

  /// @brief Return the first knot of the spline
  ///
  /// @return first knot of the spline
  const SO3& KnotsFront() const { return knots_.front(); }

  /// @brief Remove first knot of the spline and increase the start time
  void KnotsPopFront() {
    start_t_ns_ += dt_ns_;
    knots_.pop_front();
  }

  /// @brief Resize containter with knots
  ///
  /// @param[in] n number of knots
  void Resize(int n) { knots_.resize(n); }

  /// @brief Return reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return reference to the knot
  SO3& GetKnot(int i) { return knots_[i]; }

  /// @brief Return const reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return const reference to the knot
  const SO3& GetKnot(int i) const { return knots_[i]; }

  /// @brief Return const reference to deque with knots
  ///
  /// @return const reference to deque with knots
  const std::deque<SO3>& knots() const { return knots_; }

  /// @brief Return time interval in nanoseconds
  ///
  /// @return time interval in nanoseconds
  int64_t getTimeIntervalNs() const { return dt_ns_; }

  /// @brief Evaluate SO(3) B-spline
  ///
  /// @param[in] time_ns time for evaluating the value of the spline in
  /// nanoseconds
  /// @param[out] J if not nullptr, return the Jacobian of the value with
  /// respect to knots
  /// @return SO(3) value of the spline
  SO3 Evaluate(int64_t time_ns, JacobianStruct* J = nullptr) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);

    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    SO3 res = knots_[s];

    Matrix3 J_helper;

    if (J != nullptr) {
      J->start_idx = s;
      J_helper.setIdentity();
    }

    for (int i = 0; i < DEG; i++) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      Vector3 delta = r01.log();
      Vector3 kdelta = delta * coeff[i + 1];

      if (J != nullptr) {
        const Matrix3 Jl_k_delta = LeftJacobianSO3(kdelta);
        const Matrix3 Jl_inv_delta = LeftJacobianInverseSO3(delta);

        J->d_val_d_knot[i] = J_helper;
        J_helper = coeff[i + 1] * res.matrix() * Jl_k_delta * Jl_inv_delta * p0.inverse().matrix();
        J->d_val_d_knot[i] -= J_helper;
      }
      res *= SO3::exp(kdelta);
    }

    if (J != nullptr) {
      J->d_val_d_knot[DEG] = J_helper;
    }

    return res;
  }

  /// @brief Evaluate rotational velocity (first time derivative) of SO(3)
  /// B-spline in the body frame

  /// First, let's note that for scalars \f$ k, \Delta k \f$ the following
  /// holds: \f$ \exp((k+\Delta k)\phi) = \exp(k\phi)\exp(\Delta k\phi), \phi
  /// \in \mathbb{R}^3\f$. This is due to the fact that rotations around the
  /// same axis are commutative.
  ///
  /// Let's take SO(3) B-spline with N=3 as an example. The evolution in time of
  /// rotation from the body frame to the world frame is described with \f[
  ///  R_{wb}(t) = R(t) = R_i \exp(k_1(t) \log(R_{i}^{-1}R_{i+1})) \exp(k_2(t)
  ///  \log(R_{i+1}^{-1}R_{i+2})), \f] where \f$ k_1, k_2 \f$ are spline
  ///  coefficients (see detailed description of \ref So3Spline). Since
  ///  expressions under logmap do not depend on time we can rename them to
  ///  constants.
  /// \f[ R(t) = R_i \exp(k_1(t) ~ d_1) \exp(k_2(t) ~ d_2). \f]
  ///
  /// With linear approximation of the spline coefficient evolution over time
  /// \f$ k_1(t_0 + \Delta t) = k_1(t_0) + k_1'(t_0)\Delta t \f$ we can write
  /// \f{align}
  ///  R(t_0 + \Delta t) &= R_i \exp( (k_1(t_0) + k_1'(t_0) \Delta t) ~ d_1)
  ///  \exp((k_2(t_0) + k_2'(t_0) \Delta t) ~ d_2)
  ///  \\ &= R_i \exp(k_1(t_0) ~ d_1) \exp(k_1'(t_0)~ d_1 \Delta t )
  ///  \exp(k_2(t_0) ~ d_2) \exp(k_2'(t_0) ~ d_2 \Delta t )
  ///  \\ &= R_i \exp(k_1(t_0) ~ d_1)
  ///  \exp(k_2(t_0) ~ d_2) \exp(R_{a}^T k_1'(t_0)~ d_1 \Delta t )
  ///  \exp(k_2'(t_0) ~ d_2 \Delta t )
  ///  \\ &= R_i \exp(k_1(t_0) ~ d_1)
  ///  \exp(k_2(t_0) ~ d_2) \exp((R_{a}^T k_1'(t_0)~ d_1 +
  ///  k_2'(t_0) ~ d_2) \Delta t )
  ///  \\ &= R(t_0) \exp((R_{a}^T k_1'(t_0)~ d_1 +
  ///  k_2'(t_0) ~ d_2) \Delta t )
  ///  \\ &= R(t_0) \exp( \omega \Delta t ),
  /// \f} where \f$ \Delta t \f$ is small, \f$ R_{a} \in SO(3) = \exp(k_2(t_0) ~
  /// d_2) \f$ and \f$ \omega \f$ is the rotational velocity in the body frame.
  /// More explicitly we have the formula for rotational velocity in the body
  /// frame \f[ \omega = R_{a}^T k_1'(t_0)~ d_1 +  k_2'(t_0) ~ d_2. \f]
  /// Derivatives of spline coefficients can be computed with \ref
  /// BaseCoeffsWithTime similar to \ref RdSpline (detailed description). With
  /// the recursive formula computations generalize to different orders of
  /// spline N.
  ///
  /// @param[in] time_ns time for evaluating velocity of the spline in
  /// nanoseconds
  /// @return rotational velocity (3x1 vector)
  Vector3 VelocityBody(int64_t time_ns) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);
    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<1>(u, &p);
    VectorN dcoeff = pow_inv_dt_[1] * this->kCumulativeBlendingMatrix * p;

    Vector3 rot_vel;
    rot_vel.setZero();

    for (int i = 0; i < DEG; i++) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      Vector3 delta = r01.log();

      rot_vel = SO3::exp(-delta * coeff[i + 1]) * rot_vel;
      rot_vel += delta * dcoeff[i + 1];
    }

    return rot_vel;
  }

  /// @brief Evaluate rotational velocity (first time derivative) of SO(3)
  /// B-spline in the body frame
  ///
  /// @param[in] time_ns time for evaluating velocity of the spline in
  /// nanoseconds
  /// @param[out] J if not nullptr, return the Jacobian of the rotational
  /// velocity in body frame with respect to knots
  /// @return rotational velocity (3x1 vector)
  Vector3 VelocityBody(int64_t time_ns, JacobianStruct* J) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);
    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<1>(u, &p);
    VectorN dcoeff = pow_inv_dt_[1] * this->kCumulativeBlendingMatrix * p;

    Vector3 delta_vec[DEG];

    Matrix3 R_tmp[DEG];
    SO3 accum;
    SO3 exp_k_delta[DEG];

    Matrix3 Jr_delta_inv[DEG];
    Matrix3 Jr_kdelta[DEG];

    for (int i = DEG - 1; i >= 0; i--) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      delta_vec[i] = r01.log();

      Jr_delta_inv[i] = RightJacobianInverseSO3(delta_vec[i]);
      Jr_delta_inv[i] *= p1.inverse().matrix();

      Vector3 k_delta = coeff[i + 1] * delta_vec[i];
      Jr_kdelta[i] = RightJacobianSO3(-k_delta);

      R_tmp[i] = accum.matrix();
      exp_k_delta[i] = Sophus::SO3d::exp(-k_delta);
      accum *= exp_k_delta[i];
    }

    Matrix3 d_vel_d_delta[DEG];

    d_vel_d_delta[0] = dcoeff[1] * R_tmp[0] * Jr_delta_inv[0];
    Vector3 rot_vel = delta_vec[0] * dcoeff[1];
    for (int i = 1; i < DEG; i++) {
      d_vel_d_delta[i] =
          R_tmp[i - 1] * SO3::hat(rot_vel) * Jr_kdelta[i] * coeff[i + 1] + R_tmp[i] * dcoeff[i + 1];
      d_vel_d_delta[i] *= Jr_delta_inv[i];

      rot_vel = exp_k_delta[i] * rot_vel + delta_vec[i] * dcoeff[i + 1];
    }

    if (J != nullptr) {
      J->start_idx = s;
      for (int i = 0; i < N; i++) {
        J->d_val_d_knot[i].setZero();
      }
      for (int i = 0; i < DEG; i++) {
        J->d_val_d_knot[i] -= d_vel_d_delta[i];
        J->d_val_d_knot[i + 1] += d_vel_d_delta[i];
      }
    }

    return rot_vel;
  }

  /// @brief Evaluate rotational acceleration (second time derivative) of SO(3)
  /// B-spline in the body frame
  ///
  /// @param[in] time_ns time for evaluating acceleration of the spline in
  /// nanoseconds
  /// @param[out] vel_body if not nullptr, return the rotational velocity in the
  /// body frame (3x1 vector) (side computation)
  /// @return rotational acceleration (3x1 vector)
  Vector3 AccelerationBody(int64_t time_ns, Vector3* vel_body = nullptr) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);
    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<1>(u, &p);
    VectorN dcoeff = pow_inv_dt_[1] * this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<2>(u, &p);
    VectorN ddcoeff = pow_inv_dt_[2] * this->kCumulativeBlendingMatrix * p;

    SO3 r_accum;

    Vector3 rot_vel;
    rot_vel.setZero();

    Vector3 rot_accel;
    rot_accel.setZero();

    for (int i = 0; i < DEG; i++) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      Vector3 delta = r01.log();

      SO3 rot = SO3::exp(-delta * coeff[i + 1]);

      rot_vel = rot * rot_vel;
      Vector3 vel_current = dcoeff[i + 1] * delta;
      rot_vel += vel_current;

      rot_accel = rot * rot_accel;
      rot_accel += ddcoeff[i + 1] * delta + rot_vel.cross(vel_current);
    }

    if (vel_body != nullptr) {
      *vel_body = rot_vel;
    }

    return rot_accel;
  }

  /// @brief Evaluate rotational acceleration (second time derivative) of SO(3)
  /// B-spline in the body frame
  ///
  /// @param[in] time_ns time for evaluating acceleration of the spline in
  /// nanoseconds
  /// @param[out] J_accel if not nullptr, return the Jacobian of the rotational
  /// acceleration in body frame with respect to knots
  /// @param[out] vel_body if not nullptr, return the rotational velocity in the
  /// body frame (3x1 vector) (side computation)
  /// @param[out] J_vel if not nullptr, return the Jacobian of the rotational
  /// velocity in the body frame (side computation)
  /// @return rotational acceleration (3x1 vector)
  Vector3 AccelerationBody(int64_t time_ns,
                           JacobianStruct* J_accel,
                           Vector3* vel_body = nullptr,
                           JacobianStruct* J_vel = nullptr) const {
    CHECK(J_accel != nullptr);

    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);
    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<1>(u, &p);
    VectorN dcoeff = pow_inv_dt_[1] * this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<2>(u, &p);
    VectorN ddcoeff = pow_inv_dt_[2] * this->kCumulativeBlendingMatrix * p;

    Vector3 delta_vec[DEG];
    Matrix3 exp_k_delta[DEG];
    Matrix3 Jr_delta_inv[DEG];
    Matrix3 Jr_kdelta[DEG];

    Vector3 rot_vel;
    rot_vel.setZero();

    Vector3 rot_accel;
    rot_accel.setZero();

    Vector3 rot_vel_arr[DEG];
    Vector3 rot_accel_arr[DEG];

    for (int i = 0; i < DEG; i++) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      delta_vec[i] = r01.log();

      Jr_delta_inv[i] = RightJacobianInverseSO3(delta_vec[i]);
      Jr_delta_inv[i] *= p1.inverse().matrix();

      Vector3 k_delta = coeff[i + 1] * delta_vec[i];
      Jr_kdelta[i] = RightJacobianSO3(-k_delta);

      exp_k_delta[i] = Sophus::SO3d::exp(-k_delta).matrix();

      rot_vel = exp_k_delta[i] * rot_vel;
      Vector3 vel_current = dcoeff[i + 1] * delta_vec[i];
      rot_vel += vel_current;

      rot_accel = exp_k_delta[i] * rot_accel;
      rot_accel += ddcoeff[i + 1] * delta_vec[i] + rot_vel.cross(vel_current);

      rot_vel_arr[i] = rot_vel;
      rot_accel_arr[i] = rot_accel;
    }

    Matrix3 d_accel_d_delta[DEG];
    Matrix3 d_vel_d_delta[DEG];

    d_vel_d_delta[DEG - 1] =
        coeff[DEG] * exp_k_delta[DEG - 1] * SO3::hat(rot_vel_arr[DEG - 2]) * Jr_kdelta[DEG - 1] +
        Matrix3::Identity() * dcoeff[DEG];

    d_accel_d_delta[DEG - 1] =
        coeff[DEG] * exp_k_delta[DEG - 1] * SO3::hat(rot_accel_arr[DEG - 2]) * Jr_kdelta[DEG - 1] +
        Matrix3::Identity() * ddcoeff[DEG] +
        dcoeff[DEG] * (SO3::hat(rot_vel_arr[DEG - 1]) -
                       SO3::hat(delta_vec[DEG - 1]) * d_vel_d_delta[DEG - 1]);

    Matrix3 pj;
    pj.setIdentity();

    Vector3 sj;
    sj.setZero();

    for (int i = DEG - 2; i >= 0; i--) {
      sj += dcoeff[i + 2] * pj * delta_vec[i + 1];
      pj *= exp_k_delta[i + 1];

      d_vel_d_delta[i] = Matrix3::Identity() * dcoeff[i + 1];
      if (i >= 1) {
        d_vel_d_delta[i] +=
            coeff[i + 1] * exp_k_delta[i] * SO3::hat(rot_vel_arr[i - 1]) * Jr_kdelta[i];
      }

      d_accel_d_delta[i] =
          Matrix3::Identity() * ddcoeff[i + 1] +
          dcoeff[i + 1] * (SO3::hat(rot_vel_arr[i]) - SO3::hat(delta_vec[i]) * d_vel_d_delta[i]);
      if (i >= 1) {
        d_accel_d_delta[i] +=
            coeff[i + 1] * exp_k_delta[i] * SO3::hat(rot_accel_arr[i - 1]) * Jr_kdelta[i];
      }

      d_vel_d_delta[i] = pj * d_vel_d_delta[i];
      d_accel_d_delta[i] = pj * d_accel_d_delta[i] - SO3::hat(sj) * d_vel_d_delta[i];
    }

    if (J_vel != nullptr) {
      J_vel->start_idx = s;
      for (int i = 0; i < N; i++) {
        J_vel->d_val_d_knot[i].setZero();
      }
      for (int i = 0; i < DEG; i++) {
        Matrix3 val = d_vel_d_delta[i] * Jr_delta_inv[i];

        J_vel->d_val_d_knot[i] -= val;
        J_vel->d_val_d_knot[i + 1] += val;
      }
    }

    if (J_accel != nullptr) {
      J_accel->start_idx = s;
      for (int i = 0; i < N; i++) {
        J_accel->d_val_d_knot[i].setZero();
      }
      for (int i = 0; i < DEG; i++) {
        Matrix3 val = d_accel_d_delta[i] * Jr_delta_inv[i];

        J_accel->d_val_d_knot[i] -= val;
        J_accel->d_val_d_knot[i + 1] += val;
      }
    }

    if (vel_body != nullptr) {
      *vel_body = rot_vel;
    }
    return rot_accel;
  }

  /// @brief Evaluate rotational jerk (third time derivative) of SO(3)
  /// B-spline in the body frame
  ///
  /// @param[in] time_ns time for evaluating jerk of the spline in
  /// nanoseconds
  /// @param[out] vel_body if not nullptr, return the rotational velocity in the
  /// body frame (3x1 vector) (side computation)
  /// @param[out] accel_body if not nullptr, return the rotational acceleration
  /// in the body frame (3x1 vector) (side computation)
  /// @return rotational jerk (3x1 vector)
  Vector3 JerkBody(int64_t time_ns,
                   Vector3* vel_body = nullptr,
                   Vector3* accel_body = nullptr) const {
    int64_t st_ns = (time_ns - start_t_ns_);

    CHECK_GE(st_ns, 0) << fmt::format(
        "st_ns={}, time_ns={}, start_t_ns={}", st_ns, time_ns, start_t_ns_);

    const int64_t s = st_ns / dt_ns_;
    const double u = double(st_ns % dt_ns_) / double(dt_ns_);

    CHECK_GE(s, 0) << "s " << s;
    CHECK_LE(size_t(s + N), knots_.size())
        << "s " << s << " N " << N << " knots.size() " << knots_.size();

    VectorN p;
    Base::template BaseCoeffsWithTime<0>(u, &p);
    VectorN coeff = this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<1>(u, &p);
    VectorN dcoeff = pow_inv_dt_[1] * this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<2>(u, &p);
    VectorN ddcoeff = pow_inv_dt_[2] * this->kCumulativeBlendingMatrix * p;

    Base::template BaseCoeffsWithTime<3>(u, &p);
    VectorN dddcoeff = pow_inv_dt_[3] * this->kCumulativeBlendingMatrix * p;

    Vector3 rot_vel;
    rot_vel.setZero();

    Vector3 rot_accel;
    rot_accel.setZero();

    Vector3 rot_jerk;
    rot_jerk.setZero();

    for (int i = 0; i < DEG; i++) {
      const SO3& p0 = knots_[s + i];
      const SO3& p1 = knots_[s + i + 1];

      SO3 r01 = p0.inverse() * p1;
      Vector3 delta = r01.log();

      SO3 rot = SO3::exp(-delta * coeff[i + 1]);

      rot_vel = rot * rot_vel;
      Vector3 vel_current = dcoeff[i + 1] * delta;
      rot_vel += vel_current;

      rot_accel = rot * rot_accel;
      Vector3 rot_vel_cross_vel_current = rot_vel.cross(vel_current);
      rot_accel += ddcoeff[i + 1] * delta + rot_vel_cross_vel_current;

      rot_jerk = rot * rot_jerk;
      rot_jerk +=
          dddcoeff[i + 1] * delta + (ddcoeff[i + 1] * rot_vel + 2 * dcoeff[i + 1] * rot_accel -
                                     dcoeff[i + 1] * rot_vel_cross_vel_current)
                                        .cross(delta);
    }

    if (vel_body != nullptr) {
      *vel_body = rot_vel;
    }
    if (accel_body != nullptr) {
      *accel_body = rot_accel;
    }

    return rot_jerk;
  }

 private:
  std::deque<SO3> knots_;             // Knots
  CONST_FIELD(int64_t, dt_ns) = 0;    // Knot interval in nanoseconds
  FIELD(int64_t, start_t_ns) = 0;     // Start time in nanoseconds
  std::array<Scalar, 4> pow_inv_dt_;  // Array with inverse powers of dt
};

/// @brief Uniform B-spline for SE(3) of order N. Internally uses an SO(3) (\ref
/// So3Spline) spline for rotation and 3D Euclidean spline (\ref RdSpline) for
/// translation (split representaion).
/// See [[arXiv:1911.08860]](https://arxiv.org/abs/1911.08860) for more details.
///
/// @param N Order of the spline.
/// @param Scalar Underlying floating type of this spline.
template <int N, typename Scalar = double>
class Se3Spline {
 public:
  static constexpr int DEG = N - 1;  // Degree of the spline.

  using MatN = Eigen::Matrix<Scalar, N, N>;
  using VecN = Eigen::Matrix<Scalar, N, 1>;
  using VecNp1 = Eigen::Matrix<Scalar, N + 1, 1>;

  using Vec3 = Eigen::Matrix<Scalar, 3, 1>;
  using Vec6 = Eigen::Matrix<Scalar, 6, 1>;
  using Vec9 = Eigen::Matrix<Scalar, 9, 1>;
  using Vec12 = Eigen::Matrix<Scalar, 12, 1>;

  using Mat3 = Eigen::Matrix<Scalar, 3, 3>;
  using Mat6 = Eigen::Matrix<Scalar, 6, 6>;

  using Mat36 = Eigen::Matrix<Scalar, 3, 6>;
  using Mat39 = Eigen::Matrix<Scalar, 3, 9>;
  using Mat312 = Eigen::Matrix<Scalar, 3, 12>;

  using Matrix3Array = std::array<Mat3, N>;
  using Matrix36Array = std::array<Mat36, N>;
  using Matrix6Array = std::array<Mat6, N>;

  using SO3 = Sophus::SO3<Scalar>;
  using SE3 = Sophus::SE3<Scalar>;

  using PosJacobianStruct = typename RdSpline<3, N, Scalar>::JacobianStruct;
  using SO3JacobianStruct = typename So3Spline<N, Scalar>::JacobianStruct;

  /// @brief Struct to store the accelerometer residual Jacobian with
  /// respect to knots
  struct AccelPosSO3JacobianStruct {
    int start_idx;
    std::array<Mat36, N> d_val_d_knot;
  };

  /// @brief Struct to store the pose Jacobian with respect to knots
  struct PosePosSO3JacobianStruct {
    int start_idx;
    std::array<Mat6, N> d_val_d_knot;
  };

  /// @brief Constructor with knot interval and start time
  ///
  /// @param[in] time_interval_ns knot time interval in nanoseconds
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  Se3Spline(int64_t time_interval_ns, int64_t start_time_ns = 0)
      : dt_ns_(time_interval_ns),
        pos_spline_(time_interval_ns, start_time_ns),
        so3_spline_(time_interval_ns, start_time_ns) {}

  /// @brief Gererate random trajectory
  ///
  /// @param[in] n number of knots to generate
  /// @param[in] static_init if true the first N knots will be the same
  /// resulting in static initial condition
  void GenerateRandomTrajectory(int n, bool static_init = false) {
    so3_spline_.GenerateRandomTrajectory(n, static_init);
    pos_spline_.GenerateRandomTrajectory(n, static_init);
  }

  /// @brief Set the knot to particular SE(3) pose
  ///
  /// @param[in] pose SE(3) pose
  /// @param[in] i index of the knot
  void SetKnot(const Sophus::SE3d& pose, int i) {
    so3_spline_.GetKnot(i) = pose.so3();
    pos_spline_.GetKnot(i) = pose.translation();
  }

  /// @brief Reset spline to have num_knots initialized at pose
  ///
  /// @param[in] pose SE(3) pose
  /// @param[in] num_knots number of knots to initialize
  void SetKnots(const Sophus::SE3d& pose, int num_knots) {
    so3_spline_.resize(num_knots);
    pos_spline_.resize(num_knots);

    for (int i = 0; i < num_knots; i++) {
      so3_spline_.GetKnot(i) = pose.so3();
      pos_spline_.GetKnot(i) = pose.translation();
    }
  }

  /// @brief Reset spline to the knots from other spline
  ///
  /// @param[in] other spline to copy knots from
  void SetKnots(const Se3Spline<N, Scalar>& other) {
    CHECK_EQ(other.dt_ns_, dt_ns_);
    CHECK_EQ(other.pos_spline_.knots().size(), other.pos_spline_.knots().size());

    int num_knots = other.pos_spline_.knots().size();

    so3_spline_.resize(num_knots);
    pos_spline_.resize(num_knots);

    for (int i = 0; i < num_knots; i++) {
      so3_spline_.GetKnot(i) = other.so3_spline_.GetKnot(i);
      pos_spline_.GetKnot(i) = other.pos_spline_.GetKnot(i);
    }
  }

  /// @brief Add knot to the end of the spline
  ///
  /// @param[in] knot knot to add
  void KnotsPushBack(const SE3& knot) {
    so3_spline_.KnotsPushBack(knot.so3());
    pos_spline_.KnotsPushBack(knot.translation());
  }

  /// @brief Remove knot from the back of the spline
  void KnotsPopBack() {
    so3_spline_.KnotsPopBack();
    pos_spline_.KnotsPopBack();
  }

  /// @brief Return the first knot of t` he spline
  ///
  /// @return first knot of the spline
  SE3 KnotsFront() const {
    SE3 res(so3_spline_.knots().front(), pos_spline_.knots().front());

    return res;
  }

  /// @brief Remove first knot of the spline and increase the start time
  void KnotsPopFront() {
    so3_spline_.KnotsPopFront();
    pos_spline_.KnotsPopFront();

    CHECK_EQ(so3_spline_.MinTimeNs(), pos_spline_.MinTimeNs());
    CHECK_EQ(so3_spline_.knots().size(), pos_spline_.knots().size());
  }

  /// @brief Return the last knot of the spline
  ///
  /// @return last knot of the spline
  SE3 KnotBack() {
    CHECK_EQ(so3_spline_.knots().size(), pos_spline_.knots().size());
    SE3 res(so3_spline_.knots().back(), pos_spline_.knots().back());
    return res;
  }

  /// @brief Return knot with index i
  ///
  /// @param i index of the knot
  /// @return knot
  SE3 GetKnot(int i) const {
    SE3 res(GetKnotSO3(i), GetKnotPos(i));
    return res;
  }

  /// @brief Return reference to the SO(3) knot with index i
  ///
  /// @param i index of the knot
  /// @return reference to the SO(3) knot
  SO3& GetKnotSO3(int i) { return so3_spline_.GetKnot(i); }

  /// @brief Return const reference to the SO(3) knot with index i
  ///
  /// @param i index of the knot
  /// @return const reference to the SO(3) knot
  const SO3& GetKnotSO3(int i) const { return so3_spline_.GetKnot(i); }

  /// @brief Return reference to the position knot with index i
  ///
  /// @param i index of the knot
  /// @return reference to the position knot
  Vec3& GetKnotPos(int i) { return pos_spline_.GetKnot(i); }

  /// @brief Return const reference to the position knot with index i
  ///
  /// @param i index of the knot
  /// @return const reference to the position knot
  const Vec3& GetKnotPos(int i) const { return pos_spline_.GetKnot(i); }

  /// @brief Apply increment to the knot
  ///
  /// The incremernt vector consists of translational and rotational parts \f$
  /// [\upsilon, \omega]^T \f$. Given the current pose of the knot \f$ R \in
  /// SO(3), p \in \mathbb{R}^3\f$ the updated pose is: \f{align}{ R' &=
  /// \exp(\omega) R
  /// \\ p' &= p + \upsilon
  /// \f}
  ///  The increment is consistent with \ref
  /// PoseState::BoxPlusForKnot.
  ///
  /// @param[in] i index of the knot
  /// @param[in] inc 6x1 increment vector
  template <typename Derived>
  void BoxPlusForKnot(int i, const Eigen::MatrixBase<Derived>& inc) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 6);
    pos_spline_.GetKnot(i) += inc.template head<3>();
    so3_spline_.GetKnot(i) = SO3::exp(inc.template tail<3>()) * so3_spline_.GetKnot(i);
  }

  /// @brief Maximum time represented by spline
  ///
  /// @return maximum time represented by spline in nanoseconds
  int64_t MaxTimeNs() const {
    CHECK_EQ(so3_spline_.MaxTimeNs(), pos_spline_.MaxTimeNs())
        << "so3_spline.MaxTimeNs() " << so3_spline_.MaxTimeNs() << " pos_spline.MaxTimeNs() "
        << pos_spline_.MaxTimeNs();
    return pos_spline_.MaxTimeNs();
  }

  /// @brief Minimum time represented by spline
  ///
  /// @return minimum time represented by spline in nanoseconds
  int64_t MinTimeNs() const {
    CHECK_EQ(so3_spline_.MinTimeNs(), pos_spline_.MinTimeNs())
        << "so3_spline.MinTimeNs() " << so3_spline_.MinTimeNs() << " pos_spline.MinTimeNs() "
        << pos_spline_.MinTimeNs();
    return pos_spline_.MinTimeNs();
  }

  /// @brief Number of knots in the spline
  int NumKnots() const { return pos_spline_.knots().size(); }

  /// @brief Linear acceleration in the world frame.
  ///
  /// @param[in] time_ns time to evaluate linear acceleration in nanoseconds
  Vec3 LinearAccelerationWorld(int64_t time_ns) const { return pos_spline_.Acceleration(time_ns); }

  /// @brief Linear velocity in the world frame.
  ///
  /// @param[in] time_ns time to evaluate linear velocity in nanoseconds
  Vec3 LinearVelocityWorld(int64_t time_ns) const { return pos_spline_.Velocity(time_ns); }

  /// @brief Rotational velocity in the body frame.
  ///
  /// @param[in] time_ns time to evaluate rotational velocity in nanoseconds
  Vec3 RotationalVelocityBody(int64_t time_ns) const { return so3_spline_.VelocityBody(time_ns); }

  /// @brief Evaluate pose.
  ///
  /// @param[in] time_ns time to evaluate pose in nanoseconds
  /// @return SE(3) pose at time_ns
  SE3 Pose(int64_t time_ns) const {
    SE3 res;

    res.so3() = so3_spline_.Evaluate(time_ns);
    res.translation() = pos_spline_.Evaluate(time_ns);

    return res;
  }

  /// @brief Evaluate pose and compute Jacobian.
  ///
  /// @param[in] time_ns time to evaluate pose in nanoseconds
  /// @param[out] J Jacobian of the pose with respect to knots
  /// @return SE(3) pose at time_ns
  Sophus::SE3d Pose(int64_t time_ns, PosePosSO3JacobianStruct* J) const {
    Sophus::SE3d res;

    typename So3Spline<N, Scalar>::JacobianStruct Jr;
    typename RdSpline<3, N, Scalar>::JacobianStruct Jp;

    res.so3() = so3_spline_.Evaluate(time_ns, &Jr);
    res.translation() = pos_spline_.Evaluate(time_ns, &Jp);

    if (J) {
      Eigen::Matrix3d RT = res.so3().inverse().matrix();

      J->start_idx = Jr.start_idx;
      for (int i = 0; i < N; i++) {
        J->d_val_d_knot[i].setZero();
        J->d_val_d_knot[i].template topLeftCorner<3, 3>() = RT * Jp.d_val_d_knot[i];
        J->d_val_d_knot[i].template bottomRightCorner<3, 3>() = RT * Jr.d_val_d_knot[i];
      }
    }

    return res;
  }

  /// @brief Evaluate pose and compute time Jacobian.
  ///
  /// @param[in] time_ns time to evaluate pose in nanoseconds
  /// @param[out] J Jacobian of the pose with time
  void d_pose_d_t(int64_t time_ns, Vec6& J) const {
    J.template head<3>() = so3_spline_.Evaluate(time_ns).inverse() * LinearVelocityWorld(time_ns);
    J.template tail<3>() = RotationalVelocityBody(time_ns);
  }

  /// @brief Evaluate gyroscope residual.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measurement gyroscope measurement
  /// @param[in] gyro_bias gyroscope calibration
  /// @return gyroscope residual
  Vec3 GyroResidual(int64_t time_ns,
                    const Vec3& measurement,
                    const CalibGyroBias<Scalar>& gyro_bias) const {
    return so3_spline_.VelocityBody(time_ns) - gyro_bias.GetCalibrated(measurement);
  }

  /// @brief Evaluate gyroscope residual and compute Jacobians.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measurement gyroscope measurement
  /// @param[in] gyro_bias gyroscope calibration
  /// @param[out] J_knots Jacobian with respect to SO(3) spline knots
  /// @param[out] J_bias Jacobian with respect to gyroscope calibration
  /// @return gyroscope residual
  Vec3 GyroResidual(int64_t time_ns,
                    const Vec3& measurement,
                    const CalibGyroBias<Scalar>& gyro_bias,
                    SO3JacobianStruct* J_knots,
                    Mat312* J_bias = nullptr) const {
    if (J_bias != nullptr) {
      J_bias->setZero();
      J_bias->template block<3, 3>(0, 0).diagonal().array() = 1.0;
      J_bias->template block<3, 3>(0, 3).diagonal().array() = -measurement[0];
      J_bias->template block<3, 3>(0, 6).diagonal().array() = -measurement[1];
      J_bias->template block<3, 3>(0, 9).diagonal().array() = -measurement[2];
    }

    return so3_spline_.VelocityBody(time_ns, J_knots) - gyro_bias.GetCalibrated(measurement);
  }

  /// @brief Evaluate accelerometer residual.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measurement accelerometer measurement
  /// @param[in] accel_bias accelerometer calibration
  /// @param[in] g gravity
  /// @return accelerometer residual
  Vec3 AccelResidual(int64_t time_ns,
                     const Eigen::Vector3d& measurement,
                     const CalibAccelBias<Scalar>& accel_bias,
                     const Eigen::Vector3d& g) const {
    Sophus::SO3d R = so3_spline_.Evaluate(time_ns);
    Eigen::Vector3d accel_world = pos_spline_.Acceleration(time_ns);

    return R.inverse() * (accel_world + g) - accel_bias.GetCalibrated(measurement);
  }

  /// @brief Evaluate accelerometer residual and Jacobians.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measurement accelerometer measurement
  /// @param[in] accel_bias accelerometer calibration
  /// @param[in] g gravity
  /// @param[out] J_knots Jacobian with respect to spline knots
  /// @param[out] J_bias Jacobian with respect to accelerometer calibration
  /// @param[out] J_g Jacobian with respect to gravity
  /// @return accelerometer residual
  Vec3 AccelResidual(int64_t time_ns,
                     const Vec3& measurement,
                     const CalibAccelBias<Scalar>& accel_bias,
                     const Vec3& g,
                     AccelPosSO3JacobianStruct* J_knots,
                     Mat39* J_bias = nullptr,
                     Mat3* J_g = nullptr) const {
    typename So3Spline<N, Scalar>::JacobianStruct Jr;
    typename RdSpline<3, N, Scalar>::JacobianStruct Jp;

    Sophus::SO3d R = so3_spline_.Evaluate(time_ns, &Jr);
    Eigen::Vector3d accel_world = pos_spline_.Acceleration(time_ns, &Jp);

    Eigen::Matrix3d RT = R.inverse().matrix();
    Eigen::Matrix3d tmp = RT * Sophus::SO3d::hat(accel_world + g);

    CHECK_EQ(Jr.start_idx, Jp.start_idx)
        << "Jr.start_idx " << Jr.start_idx << " Jp.start_idx " << Jp.start_idx;

    CHECK_EQ(so3_spline_.knots().size(), pos_spline_.knots().size())
        << "so3_spline.knots().size() " << so3_spline_.knots().size()
        << " pos_spline.knots().size() " << pos_spline_.knots().size();

    J_knots->start_idx = Jp.start_idx;
    for (int i = 0; i < N; i++) {
      J_knots->d_val_d_knot[i].template topLeftCorner<3, 3>() = RT * Jp.d_val_d_knot[i];
      J_knots->d_val_d_knot[i].template bottomRightCorner<3, 3>() = tmp * Jr.d_val_d_knot[i];
    }

    if (J_bias != nullptr) {
      J_bias->setZero();
      J_bias->template block<3, 3>(0, 0).diagonal().array() = 1.0;
      J_bias->template block<3, 3>(0, 3).diagonal().array() = -measurement[0];
      (*J_bias)(1, 6) = -measurement[1];
      (*J_bias)(2, 7) = -measurement[1];
      (*J_bias)(2, 8) = -measurement[2];
    }

    if (J_g) {
      (*J_g) = RT;
    }

    Vec3 res = RT * (accel_world + g) - accel_bias.GetCalibrated(measurement);

    return res;
  }

  /// @brief Evaluate position residual.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measured_position position measurement
  /// @param[out] Jp if not nullptr, Jacobian with respect to knos of the
  /// position spline
  /// @return position residual
  Sophus::Vector3d PositionResidual(int64_t time_ns,
                                    const Vec3& measured_position,
                                    PosJacobianStruct* Jp = nullptr) const {
    return pos_spline_.Evaluate(time_ns, Jp) - measured_position;
  }

  /// @brief Evaluate orientation residual.
  ///
  /// @param[in] time_ns time of the measurement
  /// @param[in] measured_orientation orientation measurement
  /// @param[out] Jr if not nullptr, Jacobian with respect to knos of the
  /// SO(3) spline
  /// @return orientation residual
  Sophus::Vector3d OrientationResidual(int64_t time_ns,
                                       const SO3& measured_orientation,
                                       SO3JacobianStruct* Jr = nullptr) const {
    Sophus::Vector3d res =
        (so3_spline_.Evaluate(time_ns, Jr) * measured_orientation.inverse()).log();

    if (Jr != nullptr) {
      const Eigen::Matrix3d Jrot = LeftJacobianSO3(res);
      for (int i = 0; i < N; i++) {
        Jr->d_val_d_knot[i] = Jrot * Jr->d_val_d_knot[i];
      }
    }

    return res;
  }

  /// @brief Print knots for debugging.
  void PrintKnots() const {
    for (int i = 0; i < pos_spline_.knots().size(); i++) {
      std::cout << i << ": p:" << pos_spline_.GetKnot(i).transpose()
                << " q: " << so3_spline_.GetKnot(i).unit_quaternion().coeffs().transpose()
                << std::endl;
    }
  }

  /// @brief Print position knots for debugging.
  void PrintPosKnots() const {
    for (int i = 0; i < pos_spline_.knots().size(); i++) {
      std::cout << pos_spline_.GetKnot(i).transpose() << std::endl;
    }
  }

  /// @brief Set start time for spline
  ///
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  void set_start_time_ns(int64_t s) {
    so3_spline_.set_start_time_ns(s);
    pos_spline_.set_start_time_ns(s);
  }

 private:
  FIELD(int64_t, dt_ns) = 0;           // Knot interval in nanoseconds
  RdSpline<3, N, Scalar> pos_spline_;  // Position spline
  So3Spline<N, Scalar> so3_spline_;    // Orientation spline
};
