licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

py_library(
    name = "genpy",
    srcs = glob(["src/**/*.py"]),
    imports = ["src"],
    deps = [
        "@PyYAML//:yaml",
        "@ros_genmsg//:genmsg",
    ],
)

py_binary(
    name = "genmsg_py",
    srcs = ["scripts/genmsg_py.py"],
    main = "scripts/genmsg_py.py",
    deps = [":genpy"],
)

py_binary(
    name = "gensrv_py",
    srcs = ["scripts/gensrv_py.py"],
    main = "scripts/gensrv_py.py",
    deps = [":genpy"],
)
