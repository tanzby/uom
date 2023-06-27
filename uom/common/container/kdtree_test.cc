// Copyright @2023 UOM project. All rights reserved.

#include "uom/common/container/kdtree.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace uom {

using Eigen::Vector2d;

template <typename Tree, typename... Args>
void AssertPreOrderTraversalResultAre(const Tree& kdtree, Args&&... matchers) {
  std::vector<Vector2d> traversal_points;
  kdtree.PreOrderTraversal([&traversal_points](const typename Tree::Node* node, int /*level*/) {
    traversal_points.emplace_back(node->point);
    return true;
  });
  ASSERT_THAT(traversal_points, testing::ElementsAre(std::forward<Args>(matchers)...));
}

template <typename PointType>
const typename KdTree<PointType>::Node* UnsafeFindNodeCloseTo(const KdTree<PointType>& kdtree,
                                                              const PointType& query_point) {
  using NodeType = typename KdTree<PointType>::Node;
  const NodeType* result = nullptr;
  kdtree.PreOrderTraversal([&](const NodeType* node, int /*level*/) {
    if (node->point.isApprox(query_point)) {
      result = node;
      return false;
    }
    return true;
  });
  return result;
}

std::vector<Vector2d> MockDataPoints() {
  std::vector<Vector2d> points{
      {1, 11},
      {2, 9},
      {5, 8},
      {6, 9},
      {4, 2},
      {7, 10},
      {8, 3},
      {6, 5},
      {10, 10},
      {11, 5},
      {3, 8},
      {12, 1},
      {4, 9},
      {2, 4},
  };
  return points;
}

TEST(KdTreeTest, FindBestSplitAxisAndSortByMaxSpan) {
  KdTreeParams params;
  KdTree<Eigen::Vector3d> kdtree(params);
  std::vector<Eigen::Vector3d> mutable_points{
      {0, 1, -9},
      {1, 3, 8},
      {2, 5, 1},
      {3, 7, 5},
      {4, 9, -3},
      {5, 0, 9},
  };

  auto array_span = absl::MakeSpan(mutable_points);
  const int split_axis = kdtree.FindBestSplitAxisAndSortByMaxSpan(array_span);
  ASSERT_EQ(2, split_axis);

  ASSERT_THAT(array_span.subspan(0, array_span.size() / 2),
              testing::UnorderedElementsAre(
                  Eigen::Vector3d{0, 1, -9}, Eigen::Vector3d{4, 9, -3}, Eigen::Vector3d{2, 5, 1}));
  ASSERT_THAT(array_span.subspan(array_span.size() / 2),
              testing::UnorderedElementsAre(
                  Eigen::Vector3d{3, 7, 5}, Eigen::Vector3d{1, 3, 8}, Eigen::Vector3d{5, 0, 9}));
}

TEST(KdTreeTest, BuildTree) {
  KdTreeParams params;
  KdTree<Vector2d> kdtree(params);

  // Empty points.
  {
    const std::vector<Vector2d> points{};
    kdtree.BuildTree(points);
    ASSERT_EQ(0, kdtree.Size());
  }

  // One points.
  {
    const std::vector<Vector2d> points{
        {10, 10},
    };
    kdtree.BuildTree(points);
    ASSERT_EQ(points.size(), kdtree.Size());
  }

  // Multi-points.
  {
    const std::vector<Vector2d> points{
        {1, 11},
        {2, 9},
        {5, 8},
        {4, 2},
        {8, 3},
        {10, 10},
        {11, 5},
        {12, 1},
    };
    kdtree.BuildTree(points);
    AssertPreOrderTraversalResultAre(kdtree,
                                     Vector2d(8, 3),
                                     Vector2d(2, 9),
                                     Vector2d(5, 8),
                                     Vector2d(4, 2),
                                     Vector2d(1, 11),
                                     Vector2d(11, 5),
                                     Vector2d(12, 1),
                                     Vector2d(10, 10));
  }
}

TEST(KdTreeTest, InsertAndPushUp) {
  KdTreeParams params;
  params.max_subtree_ratio = 0.7;
  params.min_tree_size_allow_rebuild = 6;
  KdTree<Vector2d> kdtree(params);
  // using NodeType = KdTree<Vector2d>::Node;

  kdtree.Insert({5, 0});  // axis = x
  AssertPreOrderTraversalResultAre(kdtree, Vector2d(5, 0));
  const auto* node = CHECK_NOTNULL(UnsafeFindNodeCloseTo(kdtree, {5, 0}));
  ASSERT_EQ(0, node->axis);
  ASSERT_EQ(1, node->size);
  ASSERT_EQ(0, node->num_invalid_nodes);
  ASSERT_DOUBLE_EQ(0.0, node->radius);

  kdtree.Insert({6, 0});  // axis = y
  AssertPreOrderTraversalResultAre(kdtree, Vector2d(5, 0), Vector2d(6, 0));
  ASSERT_EQ(2, node->size);
  ASSERT_EQ(Vector2d(6, 0), node->range.max());
  ASSERT_EQ(Vector2d(5, 0), node->range.min());
  ASSERT_DOUBLE_EQ(1.0 / 2, node->radius);

  kdtree.Insert({4, 0});  // axis = x
  AssertPreOrderTraversalResultAre(kdtree, Vector2d(5, 0), Vector2d(4, 0), Vector2d(6, 0));
  ASSERT_EQ(Vector2d(6, 0), node->range.max());
  ASSERT_EQ(Vector2d(4, 0), node->range.min());
  ASSERT_DOUBLE_EQ(1.0, node->radius);

  kdtree.Insert({7, -1});
  AssertPreOrderTraversalResultAre(
      kdtree, Vector2d(5, 0), Vector2d(4, 0), Vector2d(6, 0), Vector2d(7, -1));

  kdtree.Insert({7, 2});
  AssertPreOrderTraversalResultAre(
      kdtree, Vector2d(5, 0), Vector2d(4, 0), Vector2d(6, 0), Vector2d(7, -1), Vector2d(7, 2));

  // The insertion call should cause a rebuilding.
  kdtree.Insert({9, 3});
  AssertPreOrderTraversalResultAre(kdtree,
                                   Vector2d(7, -1),
                                   Vector2d(5, 0),
                                   Vector2d(4, 0),
                                   Vector2d(6, 0),
                                   Vector2d(9, 3),
                                   Vector2d(7, 2));
  node = CHECK_NOTNULL(UnsafeFindNodeCloseTo(kdtree, {7, -1}));
  ASSERT_EQ(0, node->axis);
  ASSERT_EQ(6, node->size);
  ASSERT_EQ(0, node->num_invalid_nodes);
  ASSERT_EQ(Vector2d(9, 3), node->range.max());
  ASSERT_EQ(Vector2d(4, -1), node->range.min());
}

TEST(KdTreeTest, LazyDeletion) {
  KdTreeParams params;
  params.max_subtree_ratio = 0.7;
  params.max_deleted_ratio = 0.5;
  params.min_tree_size_allow_rebuild = 6;
  KdTree<Vector2d> kdtree(params);

  const std::vector<Vector2d> points = MockDataPoints();
  kdtree.BuildTree(points);
  ASSERT_EQ(points.size(), kdtree.NumNodes());

  // Remove nothing.
  using RangeType = KdTree<Vector2d>::RangeType;
  kdtree.Delete(RangeType(Vector2d{100, 100}, Vector2d{200, 200}));
  ASSERT_EQ(points.size(), kdtree.NumNodes());

  // Remove {5, 8}, {6, 9}, {7, 10}, as we use lazy-deletion, the nodes are not removed.
  // └──[6, 5](0)
  //     ├──[3, 8](1)
  //     │   ├──[2, 4](1)
  //     │   │   ├──[4, 2]
  //     │   │   └──[5, 8]
  //     │   └──[2, 9](0)
  //     │       ├──[1, 11]
  //     │       └──[4, 9]
  //     └──[6, 9](1)
  //         ├──[11, 5](0)
  //         │   ├──[8, 3]
  //         │   └──[12, 1]
  //         └──[10, 10](0)
  //             ├──[7, 10]
  kdtree.Delete(RangeType(Vector2d{4.5, 7.0}, Vector2d{7.5, 10.5}));
  ASSERT_EQ(points.size(), kdtree.NumNodes());
  ASSERT_EQ(11, kdtree.Size());
  ASSERT_TRUE(CHECK_NOTNULL(UnsafeFindNodeCloseTo(kdtree, {5, 8}))->point_deleted);
  ASSERT_TRUE(CHECK_NOTNULL(UnsafeFindNodeCloseTo(kdtree, {6, 9}))->point_deleted);
  ASSERT_TRUE(CHECK_NOTNULL(UnsafeFindNodeCloseTo(kdtree, {7, 10}))->point_deleted);

  // Remove {1, 11}, {2, 9}, {4, 9}, the tree will be rebuilt.
  // └──[6, 5](0)
  //     ├──[2, 4](1)
  //     │   ├──[4, 2]
  //     │   └──[3, 8]
  //     └──[6, 9](1)
  //         ├──[11, 5](0)
  //         │   ├──[8, 3]
  //         │   └──[12, 1]
  //         └──[10, 10](0)
  //             ├──[7, 10]
  kdtree.Delete(RangeType(Vector2d{0.9, 8.9}, Vector2d{4.1, 11.1}));
  ASSERT_EQ(10, kdtree.NumNodes());
  ASSERT_EQ(8, kdtree.Size());
}

TEST(KdTreeTest, RadiusSearch) {
  KdTreeParams params;
  params.max_subtree_ratio = 0.7;
  params.max_deleted_ratio = 0.5;
  params.min_tree_size_allow_rebuild = 6;
  KdTree<Vector2d> kdtree(params);

  const std::vector<Vector2d> points = MockDataPoints();
  kdtree.BuildTree(points);
  ASSERT_EQ(points.size(), kdtree.NumNodes());

  std::vector<Vector2d> result;

  // Case 1: no any intersection.
  kdtree.RadiusSearch({8.5, 7.0}, 2.0, &result);
  ASSERT_EQ(0, result.size());

  // Case 2: sub-tree is contained.
  kdtree.RadiusSearch({3.0, 11.0}, 2.9, &result);
  ASSERT_EQ(3, result.size());

  // Case 3: part of tree is contained.
  kdtree.RadiusSearch({6.0, 9.0}, 2.1, &result);
  ASSERT_EQ(4, result.size());
}

TEST(KdTreeTest, NearestSearch) {
  KdTreeParams params;
  params.max_subtree_ratio = 0.7;
  params.max_deleted_ratio = 0.5;
  params.min_tree_size_allow_rebuild = 6;
  KdTree<Vector2d> kdtree(params);

  const std::vector<Vector2d> points = MockDataPoints();
  kdtree.BuildTree(points);
  ASSERT_EQ(points.size(), kdtree.NumNodes());

  // We will found {{6, 5}, {5, 8}, {4, 2}, {2, 4}} with these input. To verify the lazy-deletion,
  // we remove {2, 4} and {4, 2}
  const Vector2d query_point{5, 5};
  const int k = 4;
  const KdTree<Vector2d>::RangeType deleted_range(Vector2d{1.9, 1.9}, Vector2d{4.1, 4.1});
  kdtree.Delete(deleted_range);

  // Compute expected GT.
  std::vector<std::pair<double, int>> dist_index_array;
  dist_index_array.reserve(points.size());
  for (int i = 0; i < points.size(); ++i) {
    if (deleted_range.contains(points[i])) {
      continue;
    }
    dist_index_array.emplace_back((points[i] - query_point).norm(), i);
  }
  std::sort(dist_index_array.begin(), dist_index_array.end());
  std::vector<Vector2d> expected_k_points;
  std::vector<double> expected_k_distances;
  for (int i = 0; i < k; ++i) {
    expected_k_points.emplace_back(points[dist_index_array[i].second]);
    expected_k_distances.emplace_back(dist_index_array[i].first);
  }

  // Search by kdtree.
  std::vector<Vector2d> k_points;
  std::vector<double> k_distances;
  const int k_found = kdtree.NearestSearch(query_point, k, &k_points, &k_distances);
  ASSERT_EQ(k, k_found);
  ASSERT_EQ(expected_k_points, k_points);
  ASSERT_EQ(expected_k_distances, k_distances);
}

}  // namespace uom
