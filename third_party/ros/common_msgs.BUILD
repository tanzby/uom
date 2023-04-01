licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

load(
    "@com_github_mvukov_rules_ros//ros:interfaces.bzl",
    "cc_ros_interface_library",
    "py_ros_interface_library",
    "ros_interface_library",
)

ros_interface_library(
    name = "actionlib_msgs",
    srcs = glob(["actionlib_msgs/msg/*.msg"]),
    deps = ["@ros_std_msgs//:std_msgs"],
)

cc_ros_interface_library(
    name = "cc_actionlib_msgs",
    deps = [":actionlib_msgs"],
)

py_ros_interface_library(
    name = "py_actionlib_msgs",
    deps = [":actionlib_msgs"],
)

py_binary(
    name = "genaction",
    srcs = ["actionlib_msgs/scripts/genaction.py"],
)

ros_interface_library(
    name = "diagnostic_msgs",
    srcs = glob([
        "diagnostic_msgs/msg/*.msg",
        "diagnostic_msgs/srv/*.srv",
    ]),
    deps = ["@ros_std_msgs//:std_msgs"],
)

ros_interface_library(
    name = "geometry_msgs",
    srcs = glob(["geometry_msgs/msg/*.msg"]),
    deps = ["@ros_std_msgs//:std_msgs"],
)

cc_ros_interface_library(
    name = "cc_geometry_msgs",
    deps = [":geometry_msgs"],
)

ros_interface_library(
    name = "nav_msgs",
    srcs = glob([
        "nav_msgs/action/*.action",
        "nav_msgs/msg/*.msg",
        "nav_msgs/srv/*.srv",
    ]),
    deps = [
        ":actionlib_msgs",
        ":geometry_msgs",
        "@ros_std_msgs//:std_msgs",
    ],
)

cc_ros_interface_library(
    name = "cc_nav_msgs",
    deps = [":nav_msgs"],
)

py_ros_interface_library(
    name = "py_nav_msgs",
    deps = [":nav_msgs"],
)

ros_interface_library(
    name = "sensor_msgs",
    srcs = glob([
        "sensor_msgs/msg/*.msg",
        "sensor_msgs/srv/*.srv",
    ]),
    deps = [
        ":geometry_msgs",
        "@ros_std_msgs//:std_msgs",
    ],
)

cc_ros_interface_library(
    name = "cc_sensor_msgs",
    deps = [":sensor_msgs"],
)

py_ros_interface_library(
    name = "py_sensor_msgs",
    deps = [":sensor_msgs"],
)

cc_library(
    name = "cc_sensor_msgs_lib",
    hdrs = glob(["sensor_msgs/include/sensor_msgs/**/*.h"]),
    includes = ["sensor_msgs/include"],
    deps = [":cc_sensor_msgs"],
)

py_library(
    name = "py_sensor_msgs_lib",
    srcs = glob(["sensor_msgs/src/sensor_msgs/*.py"]),
    imports = ["sensor_msgs/src"],
    deps = [":py_sensor_msgs"],
)

ros_interface_library(
    name = "visualization_msgs",
    srcs = glob([
        "visualization_msgs/msg/*.msg",
    ]),
    deps = [
        ":geometry_msgs",
        "@ros_std_msgs//:std_msgs",
    ],
)

cc_ros_interface_library(
    name = "cc_visualization_msgs",
    deps = [":visualization_msgs"],
)
