# Eigen is a C++ template library for linear algebra: vectors,
# matrices, and related algorithms.

licenses([
    # Note: Eigen is an MPL2 library that includes GPL v3 and LGPL v2.1+ code.
    #       We've taken special care to not reference any restricted code.
    "reciprocal",  # MPL2
    "notice",  # Portions BSD
])

exports_files(["COPYING.MPL2"])

cc_library(
    name = "eigen",
    hdrs = glob([
        "Eigen/*",
        "Eigen/src/**/*.h",
    ]),
    defines = select({
        "@platforms//cpu:armv7": ["EIGEN_MAX_STATIC_ALIGN_BYTES=0"],
        "//conditions:default": [],
    }) + [
        "EIGEN_MPL2_ONLY",
    ],
    includes = ["."],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "unsupported",
    hdrs = glob([
        "unsupported/Eigen/*",
        "unsupported/Eigen/src/**/*.h",
        "unsupported/Eigen/CXX11/*",
        "unsupported/Eigen/CXX11/src/**/*.h",
    ]),
    defines = [
        "EIGEN_MPL2_ONLY",
    ],
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":eigen",
    ],
)
