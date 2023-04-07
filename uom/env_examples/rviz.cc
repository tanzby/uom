// Copyright @2023 UOM project. All rights reserved.

#include <random>

#include "ros/ros.h"
#include "sensor_msgs/PointCloud.h"
#include "visualization_msgs/Marker.h"

class MockSimpleMarker {
 public:
  visualization_msgs::Marker operator()() {
    visualization_msgs::Marker marker;

    marker.id = 0;
    marker.ns = "basic_shapes";
    marker.header.frame_id = "map";
    marker.header.stamp = ros::Time::now();

    // Set the marker type.  Initially this is CUBE, and cycles between that and SPHERE, ARROW, and
    // CYLINDER
    marker.type = shape_;

    // Set the marker action.  Options are ADD, DELETE, and new in ROS Indigo: 3 (DELETEALL)
    marker.action = visualization_msgs::Marker::ADD;

    // Set the pose of the marker.  This is a full 6DOF pose relative to the frame/time specified in
    // the header
    marker.pose.position.x = 0;
    marker.pose.position.y = 0;
    marker.pose.position.z = 0;
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.pose.orientation.w = 1.0;

    // Set the scale of the marker -- 1x1x1 here means 1m on a side
    marker.scale.x = 1.0;
    marker.scale.y = 1.0;
    marker.scale.z = 1.0;

    // Set the color -- be sure to set alpha to something non-zero!
    marker.color.r = 0.0f;
    marker.color.g = 1.0f;
    marker.color.b = 0.0f;
    marker.color.a = 1.0;

    switch (shape_) {
      case visualization_msgs::Marker::CUBE:
        shape_ = visualization_msgs::Marker::SPHERE;
        break;
      case visualization_msgs::Marker::SPHERE:
        shape_ = visualization_msgs::Marker::ARROW;
        break;
      case visualization_msgs::Marker::ARROW:
        shape_ = visualization_msgs::Marker::CYLINDER;
        break;
      case visualization_msgs::Marker::CYLINDER:
        shape_ = visualization_msgs::Marker::CUBE;
        break;
    }

    return marker;
  }

 private:
  uint32_t shape_ = visualization_msgs::Marker::CUBE;
};

class MockPointCloud {
 public:
  sensor_msgs::PointCloud operator()() {
    constexpr double kCenterX = 0.0;
    constexpr double kCenterY = 0.0;
    constexpr int kSizeOfCloud = 100000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> distX(kCenterX, 2.0);
    std::normal_distribution<> distY(kCenterY, 2.0);

    sensor_msgs::PointCloud point_cloud;
    point_cloud.header.frame_id = "map";
    point_cloud.points.resize(kSizeOfCloud);

    for (int i = 0; i < kSizeOfCloud; i++) {
      const double x = distX(gen);
      const double y = distY(gen);
      point_cloud.points[i].x = x;
      point_cloud.points[i].y = y;
      point_cloud.points[i].z = std::exp(-((x * x) + (y * y)) / 4.);
    }

    sensor_msgs::ChannelFloat32 depth_channel;
    depth_channel.name = "distance";
    for (int i = 0; i < kSizeOfCloud; i++) {
      depth_channel.values.push_back(point_cloud.points[i].z);
    }
    point_cloud.channels.push_back(depth_channel);

    return point_cloud;
  }
};

int main(int argc, char** argv) {
  ros::init(argc, argv, "basic_shapes");
  ros::NodeHandle n;
  ros::Rate r(10);

  ros::Publisher marker_pub = n.advertise<visualization_msgs::Marker>("visualization_marker", 1);
  ros::Publisher cloud_pub = n.advertise<sensor_msgs::PointCloud>("pointcloud", 1);

  MockPointCloud mock_cloud;
  MockSimpleMarker mock_marker;

  while (ros::ok()) {
    cloud_pub.publish(mock_cloud());
    marker_pub.publish(mock_marker());
    r.sleep();
  }
}
