licenses(["restricted"])

# We only need to use the files in the 'lib' directory
cc_library(
    name = "lz4_lib",
    srcs = [
        "lib/lz4.c",
        "lib/lz4frame.c",
        "lib/lz4frame_static.h",
        "lib/lz4hc.c",
        "lib/xxhash.c",
        "lib/xxhash.h",
    ],
    hdrs = [
        "lib/lz4.h",
        "lib/lz4frame.h",
        "lib/lz4hc.h",
    ],
    copts = [
        "-DXXH_NAMESPACE=LZ4_",
    ],
    includes = ["lib"],
    textual_hdrs = [
        "lib/lz4.c",  # required because lz4hc.c includes "lz4.c"
    ],
    visibility = ["//visibility:public"],  # Notice: According to the license, only this target can be public
)
