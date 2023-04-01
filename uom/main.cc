#include <iostream>

#include "oneapi/tbb.h"
#include "Eigen/Core"
#include "absl/container/flat_hash_map.h"
#include "fmt/ostream.h"
#include "ceres/ceres.h"

// Hack for eigen
template <typename T>
struct fmt::formatter<T, std::enable_if_t<std:: is_base_of_v<Eigen::DenseBase<T>, T>, char>>
    : ostream_formatter {};

// test ceres.
struct CostFunctor {
  template <typename T>
  bool operator()(const T* const x, T* residual) const {
    residual[0] = 10.0 - x[0];
    return true;
  }
};

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  absl::flat_hash_map<int, Eigen::Vector3d> hashtable;
  hashtable[0] = Eigen::Vector3d(0.241, 0.542, 0.552);
  hashtable[1] = Eigen::Vector3d(4.523, 1.521, 2.125);
  fmt::print("{}: {}\n{}: {}\n", 0, hashtable[0].transpose(), 1, hashtable[1].transpose());
  
  double x = 0.5;
  const double initial_x = x;
  ceres::Problem problem;
  ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
  problem.AddResidualBlock(cost_function, nullptr, &x);

  // Run the solver!
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
  std::cout << "x : " << initial_x << " -> " << x << "\n";

  // tbb lib.
  std::cout << "Hello from oneTBB "
              << TBB_VERSION_MAJOR << "."
              << TBB_VERSION_MINOR << "."
              << TBB_VERSION_PATCH
              << "!" << std::endl;

  return 0;
}
