

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

#include "geometry_msgs/PointStamped.h"
#include "ros/ros.h"
#include "sensor_msgs/PointCloud.h"

#include "uom/common/container/kdtree.h"

using Point = sensor_msgs::PointCloud::_points_type::value_type;

std::vector<Point> MockPoints() {
  constexpr double kCenterX = 0.0;
  constexpr double kCenterY = 0.0;
  constexpr int kSizeOfCloud = 100000;

  thread_local std::random_device rd;
  thread_local std::mt19937 gen(rd());
  thread_local std::normal_distribution<> distX(kCenterX, 2.0);
  thread_local std::normal_distribution<> distY(kCenterY, 2.0);

  std::vector<Point> points;
  points.resize(kSizeOfCloud);

  for (int i = 0; i < kSizeOfCloud; i++) {
    const double x = distX(gen);
    const double y = distY(gen);
    points[i].x = x;
    points[i].y = y;
    points[i].z = std::exp(-((x * x) + (y * y)) / 2.0);
  }

  return points;
}

int main(int argc, char** argv) {
  ros::init(argc, argv, "ikdtree_demo");
  ros::NodeHandle n;
  ros::Rate r(10);

  ros::Publisher origin_cloud_pub = n.advertise<sensor_msgs::PointCloud>("origin", 1);
  ros::Publisher k_nearest_cloud_pub = n.advertise<sensor_msgs::PointCloud>("k_nearest", 1);

  sensor_msgs::PointCloud cloud_msg;
  sensor_msgs::PointCloud nearest_cloud_msg;
  cloud_msg.header.frame_id = "map";
  nearest_cloud_msg.header.frame_id = "map";
  cloud_msg.points = MockPoints();
 
  uom::KdTreeParams params;
  uom::KdTree<Eigen::Vector3d> kdtree(params);

  std::vector<Eigen::Vector3d> eigen_points;
  std::transform(cloud_msg.points.begin(),
                 cloud_msg.points.end(),
                 std::back_inserter(eigen_points),
                 [](const Point& point) { return Eigen::Vector3d(point.x, point.y, point.z); });
  kdtree.BuildTree(eigen_points);

  ros::Subscriber click_point_sub = n.subscribe<geometry_msgs::PointStamped>(
      "/clicked_point",
      1,
      [&](const boost::shared_ptr<const geometry_msgs::PointStamped> clicked_point) {
        Point query_point;
        query_point.x = clicked_point->point.x;
        query_point.y = clicked_point->point.y;
        query_point.z = 0;

        LOG(ERROR) << "click at " << query_point.x << ", " << query_point.y;

        std::vector<Eigen::Vector3d> k_points;
        kdtree.RadiusSearch(
            Eigen::Vector3d(query_point.x, query_point.y, query_point.z), 1.0, &k_points);
        nearest_cloud_msg.points.clear();
        std::transform(k_points.begin(),
                       k_points.end(),
                       std::back_inserter(nearest_cloud_msg.points),
                       [](const Eigen::Vector3d& point) {
                         Point msg_point;
                         msg_point.x = point[0];
                         msg_point.y = point[1];
                         msg_point.z = point[2];
                         return msg_point;
                       });
      });

  while (ros::ok()) {
    origin_cloud_pub.publish(cloud_msg);
    k_nearest_cloud_pub.publish(nearest_cloud_msg);
    ros::spinOnce();
    r.sleep();
  }
}
