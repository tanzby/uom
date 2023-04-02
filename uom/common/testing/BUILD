package(default_visibility = ["//visibility:public"])

cc_library(
    name = "eigen",
    hdrs = ["eigen.h"],
    deps = [
        "@eigen",
        "@gtest",
    ],
)

cc_library(
    name = "jacobian",
    hdrs = ["jacobian.h"],
    deps = [
        "@eigen",
        "@gtest",
    ],
)

cc_library(
    name = "test_main",
    srcs = ["test_main.cc"],
    deps = [
        "@gflags",
        "@glog",
        "@gtest",
    ],
)
