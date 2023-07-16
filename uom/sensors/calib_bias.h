// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include "Eigen/Dense"

#include "uom/common/defines.h"

/// @brief Static calibration for accelerometer.
///
/// Calibrates axis scaling and misalignment and has 9 parameters \f$ [b_x,
/// b_y, b_z, s_1, s_2, s_3, s_4, s_5, s_6]^T \f$.
/// \f[
/// a_c = \begin{bmatrix} s_1 + 1 & 0 & 0 \\ s_2 & s_4 + 1 & 0 \\ s_3 & s_5 &
/// s_6 + 1 \\  \end{bmatrix} a_r -  \begin{bmatrix} b_x \\ b_y \\ b_z
/// \end{bmatrix}
/// \f] where  \f$ a_c \f$ is a calibrated measurement and \f$ a_r \f$ is a
/// raw measurement. When all elements are zero applying calibration results in
/// Identity mapping.
template <typename Scalar>
class CalibAccelBias {
 public:
  using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
  using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
  using Vector9 = Eigen::Matrix<Scalar, 9, 1>;

  /// @brief  Set calibration to random values (used in unit-tests).
  void SetRandom() {
    params_.setRandom();
    params_.template head<3>() /= 10;
    params_.template tail<6>() /= 100;
  }

  /// @brief Increment the calibration vector
  ///
  /// @param inc increment vector
  void operator+=(const Vector9& inc) { params_ += inc; }

  /// @brief Return bias vector and scale matrix. See detailed description in
  /// \ref CalibAccelBias.
  void GetBiasAndScale(Vector3* accel_bias, Matrix3* accel_scale) const {
    *accel_bias = params_.template head<3>();

    accel_scale->setZero();
    accel_scale->col(0) = params_.template segment<3>(3);
    (*accel_scale)(1, 1) = params_(6);
    (*accel_scale)(2, 1) = params_(7);
    (*accel_scale)(2, 2) = params_(8);
  }

  /// @brief Calibrate the measurement. See detailed description in
  /// \ref CalibAccelBias.
  ///
  /// @param raw_measurement
  /// @return calibrated measurement
  Vector3 GetCalibrated(const Vector3& raw_measurement) const {
    Vector3 accel_bias;
    Matrix3 accel_scale;
    GetBiasAndScale(&accel_bias, &accel_scale);
    return (raw_measurement + accel_scale * raw_measurement - accel_bias);
  }

  /// @brief Invert calibration (used in unit-tests).
  ///
  /// @param calibrated_measurement
  /// @return raw measurement
  Vector3 InvertCalibration(const Vector3& calibrated_measurement) const {
    Vector3 accel_bias;
    Matrix3 accel_scale;
    GetBiasAndScale(&accel_bias, &accel_scale);
    Matrix3 accel_scale_inv = (Matrix3::Identity() + accel_scale).inverse();
    return accel_scale_inv * (calibrated_measurement + accel_bias);
  }

 private:
  FIELD(Vector9, params) = Vector9::Zero();
};

/// @brief Static calibration for gyroscope.
///
/// Calibrates rotation, axis scaling and misalignment and has 12 parameters \f$
/// [b_x, b_y, b_z, s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8, s_9]^T \f$. \f[
/// \omega_c = \begin{bmatrix} s_1 + 1 & s_4 & s_7 \\ s_2 & s_5 + 1 & s_8 \\ s_3
/// & s_6 & s_9 +1 \\  \end{bmatrix} \omega_r -  \begin{bmatrix} b_x \\ b_y
/// \\ b_z \end{bmatrix} \f] where  \f$ \omega_c \f$ is a calibrated measurement
/// and \f$ \omega_r \f$ is a raw measurement. When all elements are zero
/// applying calibration results in Identity mapping.
template <typename Scalar>
class CalibGyroBias {
 public:
  using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
  using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
  using Vector12 = Eigen::Matrix<Scalar, 12, 1>;

  /// @brief Set calibration to random values (used in unit-tests).
  void SetRandom() {
    params_.setRandom();
    params_.template head<3>() /= 10;
    params_.template tail<9>() /= 100;
  }

  /// @brief Increment the calibration vector
  ///
  /// @param inc increment vector
  void operator+=(const Vector12& inc) { params_ += inc; }

  /// @brief Return bias vector and scale matrix. See detailed description in
  /// \ref CalibGyroBias.
  void GetBiasAndScale(Vector3* gyro_bias, Matrix3* gyro_scale) const {
    *gyro_bias = params_.template head<3>();
    gyro_scale->col(0) = params_.template segment<3>(3);
    gyro_scale->col(1) = params_.template segment<3>(6);
    gyro_scale->col(2) = params_.template segment<3>(9);
  }

  /// @brief Calibrate the measurement. See detailed description in
  /// \ref CalibGyroBias.
  ///
  /// @param raw_measurement
  /// @return calibrated measurement
  Vector3 GetCalibrated(const Vector3& raw_measurement) const {
    Vector3 gyro_bias;
    Matrix3 gyro_scale;
    GetBiasAndScale(&gyro_bias, &gyro_scale);
    return (raw_measurement + gyro_scale * raw_measurement - gyro_bias);
  }

  /// @brief Invert calibration (used in unit-tests).
  ///
  /// @param calibrated_measurement
  /// @return raw measurement
  Vector3 InvertCalibration(const Vector3& calibrated_measurement) const {
    Vector3 gyro_bias;
    Matrix3 gyro_scale;
    GetBiasAndScale(&gyro_bias, &gyro_scale);
    Matrix3 gyro_scale_inv = (Matrix3::Identity() + gyro_scale).inverse();
    return gyro_scale_inv * (calibrated_measurement + gyro_bias);
  }

 private:
  FIELD(Vector12, params) = Vector12::Zero();
};
