// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include "Eigen/Core"

// Computes number of combinations that include k objects out of n.
constexpr inline int ComputeBinomialCoefficient(int n, int k) {
  if (k > n) {
    return 0;
  }
  int r = 1;
  for (int d = 1; d <= k; ++d) {
    r *= n--;
    r /= d;
  }
  return r;
}

// Given N (order of the spline), compute blending matrix for uniform B-spline evaluation.
template <int N, typename Scalar = double, bool Cumulative = false>
Eigen::Matrix<Scalar, N, N> ComputeBlendingMatrix() {
  Eigen::Matrix<double, N, N> m;
  m.setZero();

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double sum = 0;

      for (int s = j; s < N; ++s) {
        sum += std::pow(-1.0, s - j) * ComputeBinomialCoefficient(N, s - j) *
               std::pow(N - s - 1.0, N - 1.0 - i);
      }
      m(j, i) = ComputeBinomialCoefficient(N - 1, N - 1 - i) * sum;
    }
  }

  if (Cumulative) {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        m.row(i) += m.row(j);
      }
    }
  }

  int64_t factorial = 1;
  for (int i = 2; i < N; ++i) {
    factorial *= i;
  }

  return (m / factorial).template cast<Scalar>();
}

// Computes base coefficient matrix for polynomials of size N.
template <int N, typename Scalar = double>
Eigen::Matrix<Scalar, N, N> ComputeBaseCoefficients() {
  Eigen::Matrix<double, N, N> base_coefficients;

  base_coefficients.setZero();
  base_coefficients.row(0).setOnes();

  constexpr int kDeg = N - 1;
  int order = kDeg;
  for (int n = 1; n < N; n++) {
    for (int i = kDeg - order; i < N; i++) {
      base_coefficients(n, i) = (order - kDeg + i) * base_coefficients(n - 1, i);
    }
    order--;
  }
  return base_coefficients.template cast<Scalar>();
}
