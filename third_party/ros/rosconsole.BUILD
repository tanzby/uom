licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

load("@com_github_mvukov_rules_ros//ros:cc_defs.bzl", "cc_ros_library")

cc_ros_library(
    name = "rosconsole",
    srcs = [
        "src/rosconsole/impl/rosconsole_print.cpp",
        "src/rosconsole/rosconsole.cpp",
        "src/rosconsole/rosconsole_backend.cpp",
    ],
    hdrs = glob(["include/**/*.h"]),
    copts = ["-Wno-unused-parameter"],
    includes = ["include"],
    deps = [
        "@boost//:regex",
        "@boost//:system",
        "@boost//:thread",
        "@roscpp_core//:cpp_common",
        "@roscpp_core//:rostime",
    ],
)
