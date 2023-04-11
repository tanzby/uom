licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

load("@bazel_skylib//lib:dicts.bzl", "dicts")
load(
    "@com_github_mvukov_rules_ros//ros:cc_defs.bzl",
    "cc_ros_binary",
    "cc_ros_library",
)
load(
    "@com_github_mvukov_rules_ros//ros:interfaces.bzl",
    "cc_ros_interface_library",
    "ros_interface_library",
)
load(
    "@com_github_mvukov_rules_ros//third_party:expand_template.bzl",
    "expand_template",
)

cc_library(
    name = "libb64",
    srcs = glob(["utilities/xmlrpcpp/libb64/src/*.c"]),
    hdrs = glob(["utilities/xmlrpcpp/libb64/include/b64/*.h"]),
    includes = ["utilities/xmlrpcpp/libb64/include"],
    visibility = ["//visibility:private"],
)

cc_ros_library(
    name = "xmlrpcpp",
    srcs = glob(["utilities/xmlrpcpp/src/*.cpp"]),
    hdrs = glob(["utilities/xmlrpcpp/include/xmlrpcpp/*.h"]),
    includes = ["utilities/xmlrpcpp/include"],
    linkopts = ["-lm"],
    visibility = ["//visibility:public"],
    deps = [
        ":libb64",
        "@boost//:thread",
        "@roscpp_core//:cpp_common",
        "@roscpp_core//:rostime",
    ],
)

ros_interface_library(
    name = "roscpp",
    srcs = glob([
        "clients/roscpp/msg/*.msg",
        "clients/roscpp/srv/*.srv",
    ]),
)

cc_ros_interface_library(
    name = "cc_roscpp",
    deps = [":roscpp"],
)

ROS_VERSION_MAJOR = 1

ROS_VERSION_MINOR = 16

ROS_VERSION_PATCH = 0

_ROS_COMMON_H = "ros/common.h"

expand_template(
    name = "ros_common_h",
    out = _ROS_COMMON_H,
    substitutions = {
        "@roscpp_VERSION_MAJOR@": str(ROS_VERSION_MAJOR),
        "@roscpp_VERSION_MINOR@": str(ROS_VERSION_MINOR),
        "@roscpp_VERSION_PATCH@": str(ROS_VERSION_PATCH),
    },
    template = "clients/roscpp/include/ros/common.h.in",
)

_CONFIG_H = "config.h"

_CONFIG_COMMON_SUBSTITUTIONS = {
    "#cmakedefine HAVE_TRUNC": "#define HAVE_TRUNC",
    "#cmakedefine HAVE_IFADDRS_H": "#define HAVE_IFADDRS_H",
}

expand_template(
    name = "config_h",
    out = _CONFIG_H,
    substitutions = select(
        {
            "@platforms//os:linux": dicts.add(
                _CONFIG_COMMON_SUBSTITUTIONS,
                {"#cmakedefine HAVE_EPOLL": "#define HAVE_EPOLL"},
            ),
            "@platforms//os:macos": dicts.add(
                _CONFIG_COMMON_SUBSTITUTIONS,
                {"#cmakedefine HAVE_EPOLL": "/*#cmakedefine HAVE_EPOLL*/"},
            ),
        },
        no_match_error = "Only Linux and macOS are supported!",
    ),
    template = "clients/roscpp/src/libros/config.h.in",
)

cc_ros_library(
    name = "roscpp_lib",
    srcs = glob(["clients/roscpp/src/libros/**/*.cpp"]) + [_CONFIG_H],
    hdrs = glob(["clients/roscpp/include/**/*.h"]) + [_ROS_COMMON_H],
    copts = ["-Wno-unused-parameter"],
    defines = [
        "BOOST_ALLOW_DEPRECATED_HEADERS",
        "BOOST_BIND_GLOBAL_PLACEHOLDERS",
    ] + select({
        "@platforms//os:macos": ["BOOST_THREAD_HAS_CONDATTR_SET_CLOCK_MONOTONIC"],
        "//conditions:default": [],
    }),
    includes = ["clients/roscpp/include"],
    linkopts = ["-lm"],
    ros_package_name = "roscpp",
    deps = [
        ":cc_roscpp",
        ":xmlrpcpp",
        "@boost//:chrono",
        "@boost//:filesystem",
        "@boost//:scope_exit",
        "@boost//:signals2",
        "@boost//:system",
        "@ros_comm_msgs//:cc_rosgraph_msgs",
        "@ros_std_msgs//:cc_std_msgs",
        "@rosconsole",
        "@roscpp_core//:roscpp_serialization",
        "@roscpp_core//:roscpp_traits",
        "@roscpp_core//:rostime",
    ],
)

cc_ros_library(
    name = "topic_tools",
    srcs = [
        "tools/topic_tools/src/parse.cpp",
        "tools/topic_tools/src/shape_shifter.cpp",
    ],
    hdrs = glob(["tools/topic_tools/include/topic_tools/*.h"]),
    includes = ["tools/topic_tools/include"],
    deps = [
        ":roscpp_lib",
        ":xmlrpcpp",
        "@ros_std_msgs//:cc_std_msgs",
        "@rosconsole",
        "@roscpp_core//:cpp_common",
        "@roscpp_core//:roscpp_serialization",
        "@roscpp_core//:roscpp_traits",
        "@roscpp_core//:rostime",
    ],
)

cc_ros_binary(
    name = "rosout",
    srcs = ["tools/rosout/rosout.cpp"],
    deps = [
        ":roscpp_lib",
    ],
)

cc_ros_library(
    name = "message_filters",
    srcs = glob(["utilities/message_filters/src/*.cpp"]),
    hdrs = glob(["utilities/message_filters/include/**/*.h"]),
    includes = ["utilities/message_filters/include"],
    deps = [
        ":roscpp_lib",
        "@boost//:thread",
        "@rosconsole",
    ],
)

cc_ros_library(
    name = "roslz4",
    srcs = [
        "utilities/roslz4/src/lz4s.c",
        "utilities/roslz4/src/xxhash.c",
        "utilities/roslz4/src/xxhash.h",
    ],
    hdrs = glob(["utilities/roslz4/include/roslz4/*.h"]),
    copts = [
        "-DXXH_NAMESPACE=ROSLZ4_",
    ],
    includes = ["utilities/roslz4/include"],
    deps = [
        "@lz4//:lz4_lib",
        "@roscpp_core//:cpp_common",
    ],
)

cc_ros_library(
    name = "rosbag_storage",
    srcs = [
        "tools/rosbag_storage/src/bag.cpp",
        "tools/rosbag_storage/src/bag_player.cpp",
        "tools/rosbag_storage/src/buffer.cpp",
        "tools/rosbag_storage/src/bz2_stream.cpp",
        "tools/rosbag_storage/src/chunked_file.cpp",
        "tools/rosbag_storage/src/lz4_stream.cpp",
        "tools/rosbag_storage/src/message_instance.cpp",
        "tools/rosbag_storage/src/no_encryptor.cpp",
        "tools/rosbag_storage/src/query.cpp",
        "tools/rosbag_storage/src/stream.cpp",
        "tools/rosbag_storage/src/uncompressed_stream.cpp",
        "tools/rosbag_storage/src/view.cpp",
    ],
    hdrs = glob(["tools/rosbag_storage/include/rosbag/*.h"]),
    copts = ["-Wno-unused-parameter"],
    includes = ["tools/rosbag_storage/include"],
    local_defines = ["BOOST_BIND_GLOBAL_PLACEHOLDERS"],
    deps = [
        ":roslz4",
        "@boost//:filesystem",
        "@boost//:format",
        "@bzip2",
        "@console_bridge",
        "@roscpp_core//:roscpp_serialization",
        "@roscpp_core//:roscpp_traits",
        "@roscpp_core//:rostime",
    ],
)

cc_ros_library(
    name = "rosbag",
    srcs = [
        "tools/rosbag/src/player.cpp",
        "tools/rosbag/src/recorder.cpp",
        "tools/rosbag/src/time_translator.cpp",
    ],
    hdrs = glob(["tools/rosbag/include/rosbag/*.h"]),
    copts = ["-D_FILE_OFFSET_BITS=64"],
    includes = ["tools/rosbag/include"],
    deps = [
        ":rosbag_storage",
        ":roscpp_lib",
        ":topic_tools",
        "@boost//:date_time",
        "@boost//:filesystem",
        "@boost//:regex",
        "@boost//:thread",
        "@ros_comm_msgs//:cc_std_srvs",
    ],
)
