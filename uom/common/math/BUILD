package(default_visibility = ["//visibility:public"])

cc_library(
    name = "sophus_utils",
    hdrs = ["sophus_utils.h"],
    deps = [
        "@sophus",
    ],
)

cc_test(
    name = "sophus_utils_test",
    srcs = ["sophus_utils_test.cc"],
    deps = [
        ":sophus_utils",
        "//uom/common/testing:eigen",
        "//uom/common/testing:jacobian",
        "//uom/common/testing:test_main",
    ],
)

cc_library(
    name = "spline",
    hdrs = [
        "spline.h",
        "spline_common.h",
    ],
    deps = [
        "//uom/common/math:sophus_utils",
        "//uom/sensors:calib_bias",
        "@fmt",
    ],
)

cc_test(
    name = "spline_test",
    srcs = ["spline_test.cc"],
    deps = [
        ":spline",
        "//uom/common/testing:eigen",
        "//uom/common/testing:jacobian",
        "//uom/common/testing:test_main",
    ],
)
