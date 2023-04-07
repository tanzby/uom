package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "rviz",
    srcs = ["rviz.cc"],
    copts = [
        "-Wno-deprecated-declarations",
    ],
    deps = [
        "@ros_comm//:roscpp_lib",
        "@ros_common_msgs//:cc_sensor_msgs",
        "@ros_common_msgs//:cc_visualization_msgs",
    ],
)

cc_binary(
    name = "rosbag_reader",
    srcs = ["rosbag_reader.cc"],
    deps = [
        "@fmt",
        "@gflags",
        "@glog",
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

sh_test(
    name = "llvm_tools",
    srcs = ["llvm_tools.sh"],
    data = [
        "@llvm_toolchain_llvm//:bin/clang-format",
        "@llvm_toolchain_llvm//:bin/clangd",
    ],
)
