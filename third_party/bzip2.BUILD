licenses(["notice"])  # BSD

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "bzip2",
    srcs = [
        "blocksort.c",
        "bzlib.c",
        "bzlib_private.h",
        "compress.c",
        "crctable.c",
        "decompress.c",
        "huffman.c",
        "randtable.c",
    ],
    hdrs = ["bzlib.h"],
    copts = ["-Wno-unused-parameter"],
    includes = ["."],
)
