// Copyright @2023 UOM project. All rights reserved.

#include <vector>

#include "uom/common/container/kdtree.h"

#include "benchmark/benchmark.h"

namespace {

static void BM_KdTreeBuildTree(benchmark::State& state) {
  const int num_points = state.range(0);
  std::vector<Eigen::Vector3d> points;
  points.reserve(num_points);
  for (int i = 0; i < num_points; ++i) {
    points.emplace_back(50.0 * Eigen::Vector3d::Random());
  }

  uom::KdTreeParams params;
  params.min_tree_size_allow_rebuild = 10;
  params.max_deleted_ratio = 0.5;
  params.max_subtree_ratio = 2.0 / 3.0;
  uom::KdTree<Eigen::Vector3d, false> kdtree(params);
  for (auto _ : state) {
    kdtree.BuildTree(points);
    kdtree.Clear();
  }
}

static void BM_KdTreeWithObjectPoolBuildTree(benchmark::State& state) {
  const int num_points = state.range(0);
  std::vector<Eigen::Vector3d> points;
  points.reserve(num_points);
  for (int i = 0; i < num_points; ++i) {
    points.emplace_back(50.0 * Eigen::Vector3d::Random());
  }

  uom::KdTreeParams params;
  params.min_tree_size_allow_rebuild = 10;
  params.max_deleted_ratio = 0.5;
  params.max_subtree_ratio = 2.0 / 3.0;
  uom::KdTree<Eigen::Vector3d, true> kdtree(params);
  for (auto _ : state) {
    kdtree.BuildTree(points);
    kdtree.Clear();
  }
}

}  // namespace

BENCHMARK(BM_KdTreeBuildTree)->Arg(1e5);
BENCHMARK(BM_KdTreeWithObjectPoolBuildTree)->Arg(1e5);

BENCHMARK_MAIN();

// Run on (6 X 4600 MHz CPU s)
// CPU Caches:
//   L1 Data 32 KiB (x6)
//   L1 Instruction 32 KiB (x6)
//   L2 Unified 256 KiB (x6)
//   L3 Unified 9216 KiB (x1)
// Load Average: 1.05, 0.82, 0.79
// --------------------------------------------------------------------
// Benchmark                          Time             CPU   Iterations
// --------------------------------------------------------------------
// BM_KdTreeBuildTree/100000   22303385 ns     22299916 ns           31
