package(default_visibility = ["//visibility:public"])

cc_library(
    name = "kdtree",
    hdrs = ["kdtree.h"],
    deps = [
        ":object_pool",
        ":priority_queue",
        "//uom/common:defines",
        "//uom/common/utils:fmt_eigen",
        "@com_google_absl//absl/types:span",
        "@glog",
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
    name = "kdtree_benchmark",
    srcs = ["kdtree_benchmark.cc"],
    deps = [
        ":kdtree",
        "@benchmark",
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

cc_library(
    name = "priority_queue",
    hdrs = ["priority_queue.h"],
)

cc_test(
    name = "priority_queue_test",
    srcs = ["priority_queue_test.cc"],
    deps = [
        ":priority_queue",
        "//uom/common/testing:test_main",
    ],
)

cc_library(
    name = "object_pool",
    hdrs = ["object_pool.h"],
    deps = [
        "//uom/common:defines",
    ]
)

cc_test(
    name = "object_pool_test",
    srcs = ["object_pool_test.cc"],
    deps = [
        ":object_pool",
        "//uom/common/testing:test_main",
    ],
)
