licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

load(
    "@com_github_mvukov_rules_ros//ros:interfaces.bzl",
    "cc_ros_interface_library",
    "py_ros_interface_library",
    "ros_interface_library",
)

ros_interface_library(
    name = "std_msgs",
    srcs = glob(["msg/*.msg"]),
)

cc_ros_interface_library(
    name = "cc_std_msgs",
    deps = [":std_msgs"],
)

cc_library(
    name = "cc_std_msgs_headers",
    hdrs = glob(["include/std_msgs/*.h"]),
    includes = ["include"],
    deps = ["@roscpp_core"],
)

py_ros_interface_library(
    name = "py_std_msgs",
    deps = [":std_msgs"],
)
