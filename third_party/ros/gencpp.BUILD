licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

py_library(
    name = "gencpp_lib",
    srcs = glob(["src/**/*.py"]),
    imports = ["src"],
    deps = ["@ros_genmsg//:genmsg"],
)

filegroup(
    name = "templates",
    srcs = glob(["scripts/*.template"]),
)

py_binary(
    name = "gencpp",
    srcs = ["scripts/gen_cpp.py"],
    data = [":templates"],
    main = "scripts/gen_cpp.py",
    deps = [":gencpp_lib"],
)
