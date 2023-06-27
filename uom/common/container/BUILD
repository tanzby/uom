package(default_visibility = ["//visibility:public"])

cc_library(
    name = "kdtree",
    hdrs = ["kdtree.h"],
    deps = [
        "//uom/common:defines",
        "//uom/common/utils:fmt_eigen",
        "@com_google_absl//absl/types:span",
        "@gtest",
    ],
)

cc_test(
    name = "kdtree_test",
    srcs = ["kdtree_test.cc"],
    deps = [
        ":kdtree",
        "//uom/common/testing:test_main",
    ],
)

cc_binary(
    name = "kdtree_demo",
    srcs = ["kdtree_demo.cc"],
    deps = [
        ":kdtree",
        "@glog",
        "@ros_comm//:roscpp_lib",
        "@ros_common_msgs//:cc_geometry_msgs",
        "@ros_common_msgs//:cc_sensor_msgs",
    ],
)
