package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [
        "@eigen",
        "@fmt",
        "@com_google_absl//absl/container:flat_hash_map",
        "@ceres_solver//:ceres",
        "@tbb",
    ],
)
