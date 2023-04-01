licenses(["notice"])  # BSD

load("@com_github_mvukov_rules_ros//ros:cc_defs.bzl", "cc_ros_library")

package(default_visibility = ["//visibility:public"])

cc_ros_library(
    name = "cpp_common",
    srcs = glob(["cpp_common/src/*.cpp"]),
    hdrs = glob(["cpp_common/include/ros/*.h"]),
    includes = ["cpp_common/include"],
    local_defines = [
        "HAVE_EXECINFO_H",
        "HAVE_CXXABI_H",
        "HAVE_GLIBC_BACKTRACE",
    ],
    deps = [
        "@boost//:system",
        "@boost//:thread",
        "@console_bridge",
    ],
)

cc_ros_library(
    name = "rostime",
    srcs = glob(["rostime/src/*.cpp"]),
    hdrs = glob(["rostime/include/ros/**/*.h"]),
    includes = ["rostime/include"],
    deps = [
        ":cpp_common",
        "@boost//:date_time",
        "@boost//:system",
        "@boost//:thread",
    ],
)

cc_library(
    name = "roscpp_traits",
    hdrs = glob(["roscpp_traits/include/ros/*.h"]),
    includes = ["roscpp_traits/include"],
    deps = [
        ":cpp_common",
        ":rostime",
    ],
)

cc_ros_library(
    name = "roscpp_serialization",
    srcs = glob(["roscpp_serialization/src/*.cpp"]),
    hdrs = glob(["roscpp_serialization/include/ros/*.h"]),
    includes = ["roscpp_serialization/include"],
    deps = [
        ":cpp_common",
        ":roscpp_traits",
        ":rostime",
    ],
)

cc_library(
    name = "roscpp_core",
    deps = [
        ":cpp_common",
        ":roscpp_serialization",
        ":roscpp_traits",
        ":rostime",
    ],
)
