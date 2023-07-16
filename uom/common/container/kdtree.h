// Copyright @2023 UOM project. All rights reserved.

#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Eigen/Geometry"
#include "absl/types/span.h"
#include "fmt/color.h"
#include "fmt/format.h"
#include "glog/logging.h"
#include "gtest/gtest_prod.h"

#include "uom/common/container/object_pool.h"
#include "uom/common/container/priority_queue.h"
#include "uom/common/defines.h"
#include "uom/common/utils/fmt_eigen.h"

namespace uom {

struct KdTreeParams {
  // Minimal tree size allowed to rebuild.
  int min_tree_size_allow_rebuild = 10;

  // If the number of deleted nodes large then `max_deleted_ratio` * tree size of this node. The
  // tree should be rebuilt.
  double max_deleted_ratio = 0.5;

  // If the number of nodes from a sub-tree rooted a child is large then `max_subtree_ratio` * tree
  // size of its parents. The tree should be rebuilt rooted at its parent.
  double max_subtree_ratio = 2.0 / 3.0;
};

struct KdTreeSplitInfo {
  int dim = 0;
  int index = 0;
};

// `PointType` should be `Eigen::Vector<>` types. And it is not thread-safe.
template <typename PointType, bool kUseObjectPool = true>
class KdTree {
  static constexpr int kDim = PointType::RowsAtCompileTime;

 public:
  EIGEN_STATIC_ASSERT_FIXED_SIZE(PointType);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(PointType);

  using RangeType = Eigen::AlignedBox<double, kDim>;

  explicit KdTree(const KdTreeParams& params) : params_{params} {
    CHECK_GT(params.min_tree_size_allow_rebuild, 0);
  }
  ~KdTree() = default;

  // Public node type for user to access the basic field of tree node.
  struct Node;

  // Get number of valid nodes of *this.
  [[nodiscard]] int Size() const { return root_ == nullptr ? 0 : root_->size - NumInvalidNodes(); }

  // Get number of all nodes of *this.
  [[nodiscard]] int NumNodes() const { return root_ == nullptr ? 0 : root_->size; }

  // Get number of invalid nodes of *this.
  [[nodiscard]] int NumInvalidNodes() const {
    return root_ == nullptr ? 0 : root_->num_invalid_nodes;
  }

  // Erase all nodes and build the binary tree from scratch. This function is usaully called at the
  // beginning of pipeline.
  void BuildTree(const std::vector<PointType>& points);

  // Inserts a point to the tree.
  void Insert(const PointType& point);

  // Deletes points among the `range`.
  void Delete(const RangeType& range);

  // Clear the tree.
  void Clear() { root_.reset(); }

  // Radius seach from given `query_point` and `search_radius`.
  void RadiusSearch(const PointType& query_point,
                    double search_radius,
                    std::vector<PointType>* result) const;

  // Search `k` nearest neighbors of given `query_point`. It returns the number of found neighbors.
  int NearestSearch(const PointType& query_point,
                    int k,
                    std::vector<PointType>* points,
                    std::vector<double>* distances);

  // Prints the binary tree.
  void PrintTree() const;

  // The `Callback` should be `bool(const Node* node, int level)`. If returned boolean is False, the
  // traversal will be terminated.
  template <typename Callback>
  void PreOrderTraversal(Callback&& callback) const {
    PreOrderTraversalRecursive(0, root_.get(), callback);
  }

 private:
  // Node structure for kd-tree.
  struct TreeNode;

  // Structure for nearest search.
  struct NearestSearchQueueItem;

  // Pool type with custom acquirer.
  struct NodeAcquirer;
  using ObjectPool = ObjectPool<TreeNode, NodeAcquirer>;

  // Store node to a RAII strucure to avoid memory leak.
  template <bool kEnable>
  struct PtrAlias {};
  template <>
  struct PtrAlias<true> {
    using Ptr = typename ObjectPool::Ptr;
  };
  template <>
  struct PtrAlias<false> {
    using Ptr = std::unique_ptr<TreeNode>;
  };
  using NodeUniquePtr = typename PtrAlias<kUseObjectPool>::Ptr;

  // Build tree rooted at `root_` from points recursive.
  void BuildTreeToNodeRecursive(absl::Span<PointType> points, NodeUniquePtr* node);

  // Find the max span dimension and return the axis index and splitted index. Check fails if size
  // of points are less than 1.
  KdTreeSplitInfo FindBestSplitAxisAndSortByMaxSpan(absl::Span<PointType> points);

  // Deletes(lazy) nodes from trees.
  void DeleteRecursive(const RangeType& range, NodeUniquePtr* node);

  // Try to insert a point to the subtree rooted at `node`.
  void InsertRecursive(const PointType& point, NodeUniquePtr* node, int parent_axis);

  // Calculates the rebuild criterions, returns true if the sub-tree rooted at `node` needs
  // rebuilding.
  bool IsSubtreeNeedRebuild(const TreeNode& node);

  // Rebulds the sub-tree rooted at `node`.
  void Rebuild(NodeUniquePtr* node);

  // Flattens the sub-tree into the given array.
  void FlattenRecursive(const NodeUniquePtr& node, std::vector<PointType>* output) const;

  // Radius search from a node recursively.
  void RadiusSearchRecursive(const NodeUniquePtr& node,
                             const PointType& query_point,
                             double search_radius,
                             std::vector<PointType>* result) const;

  // Radius search from a node recursively.
  void NearestSearchRecursive(const NodeUniquePtr& node, const PointType& query_point, int k) const;

  // Helper functions for recursively print the tree.
  void PrintTreeRecursive(const std::string& prefix,
                          const NodeUniquePtr& node,
                          bool is_left,
                          bool is_tree_deleted) const;

  template <typename Callback>
  void PreOrderTraversalRecursive(int level, TreeNode* node, Callback&& callback) const;

  KdTreeParams params_;

  ObjectPool node_pool_;

  NodeUniquePtr root_;

  mutable PriorityQueue<NearestSearchQueueItem> nearest_search_queue_;

  FRIEND_TEST(KdTreeTest, FindBestSplitAxisAndSortByMaxSpan);

  DISALLOW_COPY_MOVE_AND_ASSIGN(KdTree);
};

template <typename PointType, bool kUseObjectPool>
struct KdTree<PointType, kUseObjectPool>::Node {
  // The point stored at *this.
  PointType point = PointType::Zero();
  // The division axis of this Node for kd-tree.
  int axis = 0;
  // The number of total nodes rooted at this node, including itself.
  int size = 1;
  // The number of total invalid nodes (being deleted).
  int num_invalid_nodes = 0;
  // The cover range of sub-tree rooted at *this.
  RangeType range;
  // The radius from center of sub-tree rooted at *this to the most edge.
  double radius = 0.0;
  // True if the sub-tree rooted at *this is deleted.
  bool tree_deleted = false;
  // True if *this is deleted.
  bool point_deleted = false;

  Node() = default;
  virtual ~Node() = default;

  // Reset all fields.
  void Reset();
};

template <typename PointType, bool kUseObjectPool>
struct KdTree<PointType, kUseObjectPool>::TreeNode : Node {
  // The left child of *this.
  NodeUniquePtr left = nullptr;
  // The right child of *this.
  NodeUniquePtr right = nullptr;
  // The parent of *this, not owned.
  TreeNode* parent = nullptr;

  TreeNode() = default;
  ~TreeNode() = default;

  // Reset all fields.
  void Reset();
  // Gether structure information from left and right sub-trees.
  void PushUp();
  // Apply structure information to left and right sub-trees.
  void PushDown();
  // Returns true if the node hasn't childs.
  bool IsLeaf() const { return left == nullptr && right == nullptr; }
  // Returns the distance from node range to given point.
  double RangeTo(const PointType& query_point) const;
};

template <typename PointType, bool kUseObjectPool>
struct KdTree<PointType, kUseObjectPool>::NodeAcquirer {
  void operator()(TreeNode* node) const { node->Reset(); }
};

template <typename PointType, bool kUseObjectPool>
struct KdTree<PointType, kUseObjectPool>::NearestSearchQueueItem {
  const TreeNode* node;
  double distance;
  bool operator<(const NearestSearchQueueItem& rhs) const { return distance >= rhs.distance; }
};

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::Node::Reset() {
  point.setZero();
  axis = 0;
  size = 1;
  num_invalid_nodes = 0;
  range.setEmpty();
  radius = 0.0;
  tree_deleted = false;
  point_deleted = false;
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::TreeNode::Reset() {
  Node::Reset();
  left = nullptr;
  right = nullptr;
  parent = nullptr;
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::TreeNode::PushUp() {
  this->size = 1;
  this->range = RangeType(this->point);
  this->num_invalid_nodes = static_cast<int>(this->point_deleted);

  auto collect_from_child = [this](TreeNode* child) {
    child->parent = this;
    this->size += child->size;
    this->range.extend(child->range);
    this->num_invalid_nodes += child->num_invalid_nodes;
  };

  if (left != nullptr) {
    collect_from_child(left.get());
  }
  if (right != nullptr) {
    collect_from_child(right.get());
  }

  this->radius = ((this->range.max() - this->range.min()) * 0.5).norm();
  this->tree_deleted =
      this->point_deleted && (!left || left->tree_deleted) && (!right || right->tree_deleted);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::TreeNode::PushDown() {
  auto apply_to_child = [this](TreeNode* child) {
    child->tree_deleted |= this->tree_deleted;
    child->point_deleted |= child->tree_deleted;
    if (this->tree_deleted) {
      child->num_invalid_nodes = child->size;
    }
  };

  if (left != nullptr) {
    apply_to_child(left.get());
  }
  if (right != nullptr) {
    apply_to_child(right.get());
  }
}

template <typename PointType, bool kUseObjectPool>
double KdTree<PointType, kUseObjectPool>::TreeNode::RangeTo(const PointType& query_point) const {
  double squared_distance = 0;
  for (int d = 0; d < kDim; ++d) {
    const double val = query_point[d];
    const double min_val = this->range.min()[d];
    const double max_val = this->range.max()[d];
    if (val <= min_val) {
      squared_distance += (val - min_val) * (val - min_val);
    } else if (val > max_val) {
      squared_distance += (val - max_val) * (val - max_val);
    }
  }
  return std::sqrt(squared_distance);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::BuildTree(const std::vector<PointType>& points) {
  if (points.empty()) {
    return;
  }
  root_.reset();

  node_pool_.Reserve(points.size());

  // Copy points as we need to sort the input points to build the tree.
  std::vector<PointType> workspace = points;
  auto workspace_span = absl::MakeSpan(workspace);
  BuildTreeToNodeRecursive(workspace_span, &root_);
}

template <typename PointType, bool kUseObjectPool>
KdTreeSplitInfo KdTree<PointType, kUseObjectPool>::FindBestSplitAxisAndSortByMaxSpan(
    absl::Span<PointType> points) {
  CHECK_GT(points.size(), 1);

  std::array<double, kDim> min_value;
  std::array<double, kDim> max_value;
  for (int d = 0; d < kDim; ++d) {
    min_value[d] = points[0][d];
    max_value[d] = points[0][d];
  }
  for (int i = 1; i < static_cast<int>(points.size()); ++i) {
    for (int d = 0; d < kDim; ++d) {
      min_value[d] = std::min(min_value[d], points[i][d]);
      max_value[d] = std::max(max_value[d], points[i][d]);
    }
  }

  KdTreeSplitInfo split;
  split.dim = 0;
  double max_span = max_value[0] - min_value[0];
  for (int d = 1; d < kDim; ++d) {
    const double curr_span = max_value[d] - min_value[d];
    if (curr_span > max_span) {
      split.dim = d;
      max_span = curr_span;
    }
  }

  split.index = points.size() / 2;
  auto cmp = [dim = split.dim](auto&& lhs, auto&& rhs) { return lhs[dim] < rhs[dim]; };
  std::nth_element(points.begin(), points.begin() + split.index, points.end(), cmp);

  // Find the element that all right's neighhors are larger but not equal to it, so that all the
  // nodes from left sub-tree are smaller that root.
  while (split.index > 0 && points[split.index - 1] == points[split.index]) {
    --split.index;
  }

  return split;
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::BuildTreeToNodeRecursive(absl::Span<PointType> points,
                                                                 NodeUniquePtr* root_node) {
  CHECK(root_node != nullptr);
  if (points.empty()) {
    return;
  }

  // Create new node if not exist.
  if (*root_node == nullptr) {
    if constexpr (kUseObjectPool) {
      *root_node = node_pool_.Acquire();
    } else {
      *root_node = std::make_unique<TreeNode>();
    }
  }
  TreeNode& root = (*root_node->get());

  // Recursively build sub-trees.
  if (points.size() > 1) {
    const KdTreeSplitInfo split = FindBestSplitAxisAndSortByMaxSpan(points);
    root.axis = split.dim;
    root.point = points[split.index];
    BuildTreeToNodeRecursive(points.subspan(0, split.index), &root.left);
    BuildTreeToNodeRecursive(points.subspan(split.index + 1), &root.right);
  } else {
    root.point = points[0];
  }

  root.PushUp();
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::Insert(const PointType& point) {
  InsertRecursive(point, &root_, root_ == nullptr ? 0 : root_->axis);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::Delete(const RangeType& range) {
  DeleteRecursive(range, &root_);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::DeleteRecursive(const RangeType& range,
                                                        NodeUniquePtr* node) {
  if (node == nullptr || (*node) == nullptr || (*node)->tree_deleted) {
    return;
  }

  TreeNode* node_ptr = (*node).get();
  node_ptr->PushDown();

  // Case 1: no any overlap range.
  if (!range.intersects(node_ptr->range)) {
    return;
  }

  // Case 2: the total sub-tree is contained by deleted range. Note that, only the `tree_deleted`
  // label is set but the `point_deleted` field of childrent is not.
  if (range.contains(node_ptr->range)) {
    node_ptr->tree_deleted = true;
    node_ptr->point_deleted = true;
    node_ptr->num_invalid_nodes = node_ptr->size;
    return;
  }

  // Case 3: the point of node is contained by deleted range.
  if (!node_ptr->point_deleted && range.contains(node_ptr->point)) {
    node_ptr->point_deleted = true;
  }
  DeleteRecursive(range, &node_ptr->left);
  DeleteRecursive(range, &node_ptr->right);
  node_ptr->PushUp();

  if (IsSubtreeNeedRebuild(*node_ptr)) {
    Rebuild(node);
  }
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::InsertRecursive(const PointType& point,
                                                        NodeUniquePtr* node,
                                                        int parent_axis) {
  if ((*node) == nullptr) {
    if constexpr (kUseObjectPool) {
      *node = node_pool_.Acquire();
    } else {
      *node = std::make_unique<TreeNode>();
    }
    (*node)->point = point;
    (*node)->axis = (node == &root_) ? 0 : ((parent_axis + 1) % kDim);
    (*node)->PushUp();
    return;
  }

  const int curr_axis = ((*node)->axis);
  if (point[curr_axis] < (*node)->point[curr_axis]) {
    InsertRecursive(point, &(*node)->left, curr_axis);
  } else {
    InsertRecursive(point, &(*node)->right, curr_axis);
  }

  (*node)->PushUp();

  if (IsSubtreeNeedRebuild((*node->get()))) {
    Rebuild(node);
  }
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::Rebuild(NodeUniquePtr* node) {
  if (node == nullptr || (*node) == nullptr) {
    return;
  }
  std::vector<PointType> points;
  points.reserve((*node)->size);
  FlattenRecursive(*node, &points);
  node->release();
  BuildTreeToNodeRecursive(absl::MakeSpan(points), node);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::FlattenRecursive(const NodeUniquePtr& node,
                                                         std::vector<PointType>* output) const {
  if (node == nullptr || node->tree_deleted) {
    return;
  }
  if (!node->point_deleted) {
    output->emplace_back(node->point);
  }
  FlattenRecursive(node->left, output);
  FlattenRecursive(node->right, output);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::RadiusSearch(const PointType& query_point,
                                                     double search_radius,
                                                     std::vector<PointType>* result) const {
  CHECK_GE(search_radius, 0.0);
  CHECK_NOTNULL(result)->clear();
  RadiusSearchRecursive(root_, query_point, search_radius, result);
}

template <typename PointType, bool kUseObjectPool>
int KdTree<PointType, kUseObjectPool>::NearestSearch(const PointType& query_point,
                                                     int k,
                                                     std::vector<PointType>* points,
                                                     std::vector<double>* distances) {
  CHECK_GE(k, 1);
  CHECK_NOTNULL(points)->clear();
  CHECK_NOTNULL(distances)->clear();

  // Performs recursively searching.
  nearest_search_queue_.clear();
  NearestSearchRecursive(root_, query_point, k);
  const int k_found = std::min<int>(k, nearest_search_queue_.size());
  if (k_found == 0) {
    return 0;
  }

  // Collects seach results.
  points->resize(k_found);
  distances->resize(k_found);
  for (int i = 0; i < k_found; i++) {
    (*points)[k_found - i - 1] = nearest_search_queue_.top().node->point;
    (*distances)[k_found - i - 1] = nearest_search_queue_.top().distance;
    nearest_search_queue_.pop();
  }
  return k_found;
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::RadiusSearchRecursive(
    const NodeUniquePtr& node,
    const PointType& query_point,
    double search_radius,
    std::vector<PointType>* result) const {
  if (node == nullptr) {
    return;
  }

  // Force update childs' status (point_deleted).
  const_cast<TreeNode*>(node.get())->PushDown();  // NOLINT

  const PointType range_center = (node->range.min() + node->range.max()) * 0.5;
  const double tree_to_query_dist = (range_center - query_point).norm();

  // Case 1: No possible to find a node from tree in search radius.
  if (tree_to_query_dist > search_radius + node->radius) {
    return;
  }

  // Case 2: The `search radius` range contains the sub-tree rooted at `node`.
  if (tree_to_query_dist <= search_radius - node->radius) {
    FlattenRecursive(node, result);
    return;
  }

  // Case 3: The circle with `search radius` intersects the sub-tree rooted at `node`.
  if (!node->point_deleted && (node->point - query_point).norm() <= search_radius) {
    result->emplace_back(node->point);
  }
  RadiusSearchRecursive(node->left, query_point, search_radius, result);
  RadiusSearchRecursive(node->right, query_point, search_radius, result);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::NearestSearchRecursive(const NodeUniquePtr& node,
                                                               const PointType& query_point,
                                                               int k) const {
  if (node == nullptr || node->tree_deleted) {
    return;
  }

  // Helper function to get the distance of fathest neighbor in queue.
  auto max_distance = [this]() { return nearest_search_queue_.top().distance; };

  if (!node->point_deleted) {
    const double dist = (node->point - query_point).norm();
    if (nearest_search_queue_.size() < k || dist < max_distance()) {
      if (nearest_search_queue_.size() >= k) {
        nearest_search_queue_.pop();
      }
      nearest_search_queue_.emplace(NearestSearchQueueItem{node.get(), dist});
    }
  }

  const double to_left_range_distance =
      node->left ? node->left->RangeTo(query_point) : std::numeric_limits<double>::max();
  const double to_right_range_distance =
      node->right ? node->right->RangeTo(query_point) : std::numeric_limits<double>::max();

  if (nearest_search_queue_.size() < k ||
      (to_left_range_distance < max_distance() && to_right_range_distance < max_distance())) {
    // Search left first then right, or right first then left.
    if (to_left_range_distance < to_right_range_distance) {
      NearestSearchRecursive(node->left, query_point, k);
      if (nearest_search_queue_.size() > k || to_right_range_distance < max_distance()) {
        NearestSearchRecursive(node->right, query_point, k);
      }
    } else {  // to_left_range_distance >= to_right_range_distance
      NearestSearchRecursive(node->right, query_point, k);
      if (nearest_search_queue_.size() > k || to_left_range_distance < max_distance()) {
        NearestSearchRecursive(node->left, query_point, k);
      }
    }
    return;
  }

  // If the queue is full, check the range first.
  if (to_left_range_distance < max_distance()) {
    NearestSearchRecursive(node->left, query_point, k);
  }
  if (to_right_range_distance < max_distance()) {
    NearestSearchRecursive(node->right, query_point, k);
  }
}

template <typename PointType, bool kUseObjectPool>
bool KdTree<PointType, kUseObjectPool>::IsSubtreeNeedRebuild(const TreeNode& node) {
  if (node.size < params_.min_tree_size_allow_rebuild) {
    return false;
  }

  TreeNode* child_node = node.left == nullptr ? node.right.get() : node.left.get();
  const double subtree_ratio = static_cast<double>(child_node->size) / (node.size - 1);
  const double deleted_ratio = static_cast<double>(node.num_invalid_nodes) / node.size;

  if (subtree_ratio > params_.max_subtree_ratio ||
      subtree_ratio < 1.0 - params_.max_subtree_ratio) {
    return true;
  }
  return deleted_ratio > params_.max_deleted_ratio;
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::PrintTree() const {
  PrintTreeRecursive("", root_, false, false);
}

template <typename PointType, bool kUseObjectPool>
void KdTree<PointType, kUseObjectPool>::PrintTreeRecursive(const std::string& prefix,
                                                           const NodeUniquePtr& node,
                                                           bool is_left,
                                                           bool is_tree_deleted) const {
  static Eigen::IOFormat kCommaInitFmt(
      Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");

  if (node == nullptr) {
    return;
  }
  is_tree_deleted |= node->tree_deleted;
  const bool is_deleted = node->point_deleted || is_tree_deleted;
  const auto node_color = fmt::fg(is_deleted ? fmt::color::red : fmt::color::green);
  std::stringstream sstream;
  sstream << node->point.transpose().format(kCommaInitFmt);
  const std::string axis = node->IsLeaf() ? std::string() : fmt::format("({})", node->axis);
  fmt::print("{}{}{}{}\n",
             prefix,
             (is_left ? "├──" : "└──"),
             fmt::styled(sstream.str(), node_color),
             fmt::styled(axis, fmt::fg(fmt::color::gray)));
  const auto next_prefix = prefix + (is_left ? "│   " : "    ");
  PrintTreeRecursive(next_prefix, node->left, true, is_tree_deleted);
  PrintTreeRecursive(next_prefix, node->right, false, is_tree_deleted);
}

template <typename PointType, bool kUseObjectPool>
template <typename Callback>
void KdTree<PointType, kUseObjectPool>::PreOrderTraversalRecursive(int level,
                                                                   TreeNode* node,
                                                                   Callback&& callback) const {
  if (node == nullptr) {
    return;
  }
  if (!callback(static_cast<const Node*>(node), level)) {
    return;
  }
  PreOrderTraversalRecursive(level + 1, node->left.get(), callback);
  PreOrderTraversalRecursive(level + 1, node->right.get(), callback);
}

}  // namespace uom
