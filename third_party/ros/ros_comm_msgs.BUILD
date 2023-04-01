licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

load(
    "@com_github_mvukov_rules_ros//ros:interfaces.bzl",
    "cc_ros_interface_library",
    "py_ros_interface_library",
    "ros_interface_library",
)

ros_interface_library(
    name = "rosgraph_msgs",
    srcs = glob(["rosgraph_msgs/msg/*.msg"]),
    deps = ["@ros_std_msgs//:std_msgs"],
)

ros_interface_library(
    name = "std_srvs",
    srcs = glob(["std_srvs/srv/*.srv"]),
)

cc_ros_interface_library(
    name = "cc_rosgraph_msgs",
    deps = [":rosgraph_msgs"],
)

cc_ros_interface_library(
    name = "cc_std_srvs",
    deps = [":std_srvs"],
)

py_ros_interface_library(
    name = "py_rosgraph_msgs",
    deps = [":rosgraph_msgs"],
)

py_ros_interface_library(
    name = "py_std_srvs",
    deps = [":std_srvs"],
)
