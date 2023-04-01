package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "publish_rviz_marker",
    srcs = ["publish_rviz_marker.cc"],
    deps = [
        "@ros_comm//:roscpp_lib",
        "@ros_common_msgs//:cc_visualization_msgs",
    ],
)

cc_binary(
    name = "rosbag_reader",
    srcs = ["rosbag_reader.cc"],
    deps = [
        "@glog",
        "@fmt",
        "@gflags",
        "@ros_comm//:rosbag_storage",
        "@ros_common_msgs//:cc_sensor_msgs",
    ],
)

cc_binary(
    name = "sophus",
    srcs = ["sophus.cc"],
    deps = [
        "@sophus",
    ],
)
