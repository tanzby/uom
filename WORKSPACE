workspace(name = "uom")

load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Hedron's Compile Commands Extractor for Bazel
# https://github.com/hedronvision/bazel-compile-commands-extractor
http_archive(
    name = "hedron_compile_commands",
    sha256 = "7f6ebb62298694d8cf3ecaed81b3bb48de559819ac1909d4055abdc8c0ae1000",
    strip_prefix = "bazel-compile-commands-extractor-800b9cd260ce3878e94abb7d583a7c0865f7d967",
    url = "https://github.com/hedronvision/bazel-compile-commands-extractor/archive/800b9cd260ce3878e94abb7d583a7c0865f7d967.tar.gz",
)

load("@hedron_compile_commands//:workspace_setup.bzl", "hedron_compile_commands_setup")

hedron_compile_commands_setup()

new_git_repository(
    name = "eigen",
    build_file = "//third_party/eigen:eigen.BUILD",
    commit = "3147391d946bb4b6c68edd901f2add6ac1f31f8c",  # 3.4.0
    patch_args = ["-p1"],
    patches = [
        "//third_party/eigen:0001-delete-unused-but-set-variables.patch",
    ],
    remote = "https://gitlab.com/libeigen/eigen.git",
)

new_git_repository(
    name = "ceres_solver",
    build_file = "//third_party:ceres.BUILD",
    remote = "https://github.com/ceres-solver/ceres-solver.git",
    tag = "2.1.0",
)

new_git_repository(
    name = "zlib",
    build_file = "//third_party:zlib.BUILD",
    remote = "https://github.com/madler/zlib.git",
    tag = "v1.2.13",
)

git_repository(
    name = "com_google_absl",
    remote = "https://github.com/abseil/abseil-cpp.git",
    repo_mapping = {
        "@com_google_googletest": "@gtest",
        "@com_github_google_benchmark": "@benchmark",
    },
    tag = "20230125.2",
)

git_repository(
    name = "fmt",
    commit = "a33701196adfad74917046096bf5a2aa0ab0bb50",  # 9.1.0
    patch_cmds = [
        "mv support/bazel/.bazelrc .bazelrc",
        "mv support/bazel/.bazelversion .bazelversion",
        "mv support/bazel/BUILD.bazel BUILD.bazel",
        "mv support/bazel/WORKSPACE.bazel WORKSPACE.bazel",
    ],
    remote = "https://github.com/fmtlib/fmt.git",
)

git_repository(
    name = "bazel_skylib",
    commit = "141432789c92e9db2402ef0be58e2a2d2c4dd1fd",
    remote = "https://github.com/bazelbuild/bazel-skylib.git",
)

git_repository(
    name = "gflags",
    commit = "e171aa2d15ed9eb17054558e0b3a6a413bb01067",  # v2.2.2
    remote = "https://github.com/gflags/gflags.git",
)

git_repository(
    name = "glog",
    commit = "b33e3bad4c46c8a6345525fd822af355e5ef9446",  # v0.6.0
    patch_args = ["-p1"],
    patches = [
        "//third_party:glog.patch",
    ],
    remote = "https://github.com/google/glog.git",
    repo_mapping = {
        "@com_github_gflags_gflags": "@gflags",
    },
)

git_repository(
    name = "com_google_protobuf",
    remote = "https://github.com/google/protobuf.git",
    tag = "v3.19.6",
)

git_repository(
    name = "rules_python",
    commit = "c394c46fc1b21853bc68a7a47c1fe2db828d1dd0",
    remote = "https://github.com/bazelbuild/rules_python.git",
)

git_repository(
    name = "tbb",
    patch_args = ["-p1"],
    patches = ["//third_party:tbb.patch"],
    remote = "https://github.com/oneapi-src/oneTBB.git",
    tag = "v2021.8.0",
)

git_repository(
    name = "openexr",
    remote = "https://github.com/AcademySoftwareFoundation/openexr.git",
    repo_mapping = {
        "@net_zlib_zlib": "@zlib",
    },
    tag = "v3.1.6",
)

http_archive(
    name = "gtest",
    sha256 = "82ad62a4e26c199de52a707778334e80f6b195dd298d48d520d8507d2bcb88c4",
    strip_prefix = "googletest-2d4f208765af7fa376b878860a7677ecc0bc390a",
    urls = [
        "https://github.com/google/googletest/archive/2d4f208765af7fa376b878860a7677ecc0bc390a.zip",
    ],
)

http_archive(
    name = "benchmark",
    sha256 = "ede6830512f21490eeea1f238f083702eb178890820c14451c1c3d69fd375b19",
    strip_prefix = "benchmark-a3235d7b69c84e8c9ff8722a22b8ac5e1bc716a6",
    urls = [
        "https://github.com/google/benchmark/archive/a3235d7b69c84e8c9ff8722a22b8ac5e1bc716a6.zip",
    ],
)

http_archive(
    name = "sophus",
    build_file = "//third_party:sophus.BUILD",
    sha256 = "1adc0e083e0c24abe5eb92e8613f928eff074db27dda533c285b9203e1e053b3",
    strip_prefix = "Sophus-61f9a9815f7f5d4d9dcb7f4ad9f4f42ab3563108",
    urls = [
        "https://github.com/strasdat/Sophus/archive/61f9a9815f7f5d4d9dcb7f4ad9f4f42ab3563108.zip",
    ],
)

http_archive(
    name = "Imath",
    build_file = "//third_party:Imath.BUILD",
    sha256 = "bff1fa140f4af0e7f02c6cb78d41b9a7d5508e6bcdfda3a583e35460eb6d4b47",
    strip_prefix = "Imath-3.1.7",
    urls = [
        "https://github.com/AcademySoftwareFoundation/Imath/archive/refs/tags/v3.1.7.tar.gz",
    ],
)

http_archive(
    name = "nasm",
    build_file = "//third_party:nasm.BUILD",
    sha256 = "63ec86477ad3f0f6292325fd89e1d93aea2e2fd490070863f17d48f7cd387011",
    strip_prefix = "nasm-2.13.03",
    urls = [
        "http://www.nasm.us/pub/nasm/releasebuilds/2.13.03/nasm-2.13.03.tar.bz2",
    ],
)

http_archive(
    name = "libjpeg_turbo",
    build_file = "//third_party:libjpeg_turbo.BUILD",
    sha256 = "a78b05c0d8427a90eb5b4eb08af25309770c8379592bb0b8a863373128e6143f",
    strip_prefix = "libjpeg-turbo-2.1.4",
    urls = [
        "https://github.com/libjpeg-turbo/libjpeg-turbo/archive/refs/tags/2.1.4.tar.gz",
    ],
)

http_archive(
    name = "com_github_nelhage_rules_boost",
    repo_mapping = {
        "@org_lzma_lzma": "@xz",
        "@com_github_facebook_zstd": "@zstd",
        "@org_bzip_bzip2": "@bzip2",
        "@net_zlib_zlib": "@zlib",
    },
    sha256 = "d94689a734828e5cc9f6624e58c37cb201b87c3068669bc3d2d0d6400c421667",
    strip_prefix = "rules_boost-db7d2da158e2dc49c5ba8007da4b4e3dd087b57a",
    url = "https://github.com/nelhage/rules_boost/archive/db7d2da158e2dc49c5ba8007da4b4e3dd087b57a.tar.gz",
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()

http_archive(
    name = "lz4",
    build_file = "//third_party:lz4.BUILD",
    sha256 = "030644df4611007ff7dc962d981f390361e6c97a34e5cbc393ddfbe019ffe2c1",
    strip_prefix = "lz4-1.9.3",
    urls = [
        "https://github.com/lz4/lz4/archive/v1.9.3.tar.gz",
    ],
)

http_archive(
    name = "zstd",
    build_file = "//third_party:zstd.BUILD",
    sha256 = "98e9c3d949d1b924e28e01eccb7deed865eefebf25c2f21c702e5cd5b63b85e1",
    strip_prefix = "zstd-1.5.5",
    url = "https://github.com/facebook/zstd/archive/v1.5.5/zstd-1.5.5.tar.gz",
    urls = [
        "https://github.com/facebook/zstd/archive/v1.5.5.tar.gz",
    ],
)

http_archive(
    name = "xz",
    build_file = "//third_party:xz.BUILD",
    sha256 = "0d2b89629f13dd1a0602810529327195eff5f62a0142ccd65b903bc16a4ac78a",
    strip_prefix = "xz-5.2.5",
    urls = [
        "https://github.com/xz-mirror/xz/archive/v5.2.5.tar.gz",
    ],
)

http_archive(
    name = "libwebp",
    build_file = "//third_party:libwebp.BUILD",
    sha256 = "01bcde6a40a602294994050b81df379d71c40b7e39c819c024d079b3c56307f4",
    strip_prefix = "libwebp-1.2.1",
    urls = [
        "https://storage.googleapis.com/mirror.tensorflow.org/github.com/webmproject/libwebp/archive/v1.2.1.tar.gz",
        "https://github.com/webmproject/libwebp/archive/v1.2.1.tar.gz",
    ],
)

http_archive(
    name = "libjpeg",
    build_file = "//third_party:libjpeg.BUILD",
    sha256 = "566241ad815df935390b341a5d3d15a73a4000e5aab40c58505324c2855cbbb8",
    strip_prefix = "jpeg-9b",
    urls = [
        "http://www.ijg.org/files/jpegsrc.v9b.tar.gz",
    ],
)

http_archive(
    name = "libtiff",
    build_file = "//third_party:libtiff.BUILD",
    sha256 = "c7a1d9296649233979fa3eacffef3fa024d73d05d589cb622727b5b08c423464",
    strip_prefix = "tiff-4.5.0",
    urls = [
        "https://download.osgeo.org/libtiff/tiff-4.5.0.tar.gz",
    ],
)

http_archive(
    name = "bzip2",
    build_file = "//third_party:bzip2.BUILD",
    sha256 = "ab5a03176ee106d3f0fa90e381da478ddae405918153cca248e682cd0c4a2269",
    strip_prefix = "bzip2-1.0.8",
    urls = [
        "https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz",
    ],
)

http_archive(
    name = "empy",
    build_file = "//third_party:empy.BUILD",
    sha256 = "9841e36dd26c7f69fe1005f9d9e078e41bdd50dd56fc77837ae390fb6af1aed7",
    strip_prefix = "empy-3.3.3",
    urls = [
        "https://mirror.bazel.build/www.alcyone.com/software/empy/empy-3.3.3.tar.gz",
        "http://www.alcyone.com/software/empy/empy-3.3.3.tar.gz",
    ],
)

http_archive(
    name = "pangolin",
    build_file = "//third_party:pangolin.BUILD",
    sha256 = "54aeb7c817711eb76e776c638962fb315bab4b3f57f6b91c62de3c37f0512f4b",
    strip_prefix = "Pangolin-fe57db532ba2a48319f09a4f2106cc5625ee74a9",
    urls = [
        "https://github.com/stevenlovegrove/Pangolin/archive/fe57db532ba2a48319f09a4f2106cc5625ee74a9.tar.gz",
    ],
)

new_local_repository(
    name = "opengl",
    build_file = "//third_party:opengl.BUILD",
    path = "/usr/include",
)

new_local_repository(
    name = "glew",
    build_file = "//third_party:glew.BUILD",
    path = "/usr/include",
)

http_archive(
    name = "com_github_mvukov_rules_ros",
    sha256 = "d2d96acc8702c4eea2c1560da874bc87be7f2c8dd5c751da598a065bec3410a7",
    strip_prefix = "rules_ros-7d71becdf50a771a5df0ab20934fc52b33ef605e",
    urls = [
        "https://github.com/mvukov/rules_ros/archive/7d71becdf50a771a5df0ab20934fc52b33ef605e.zip",
    ],
)

http_archive(
    name = "console_bridge",
    build_file = "//third_party/ros:console_bridge.BUILD",
    patch_args = ["-p1"],
    patches = ["//third_party/ros:console_bridge.patch"],
    sha256 = "2ff175a9bb2b1849f12a6bf972ce7e4313d543a2bbc83b60fdae7db6e0ba353f",
    strip_prefix = "console_bridge-1.0.1",
    urls = [
        "https://github.com/ros/console_bridge/archive/1.0.1.tar.gz",
    ],
)

http_archive(
    name = "roscpp_core",
    build_file = "//third_party/ros:roscpp_core.BUILD",
    sha256 = "a2aa77814ed97b48995c872a405c51f6b0f1ab9d40e38ece483852bbd273ad7b",
    strip_prefix = "roscpp_core-0.7.2",
    urls = [
        "https://github.com/ros/roscpp_core/archive/0.7.2.tar.gz",
    ],
)

http_archive(
    name = "ros_genmsg",
    build_file = "//third_party/ros:genmsg.BUILD",
    sha256 = "0e414846823a2aaa7781f81268251c7c9a45ff96cef8e6a78bbbbcf7e4c28d56",
    strip_prefix = "genmsg-0.5.16",
    urls = [
        "https://github.com/ros/genmsg/archive/0.5.16.tar.gz",
    ],
)

http_archive(
    name = "ros_gencpp",
    build_file = "//third_party/ros:gencpp.BUILD",
    patch_args = ["-p1"],
    patches = ["//third_party:ros_gencpp.patch"],
    sha256 = "05acfeeb1bbc374356bf7674fee2a7aab3bf6a48ebad4a06fd0f0d4455a60720",
    strip_prefix = "gencpp-0.6.5",
    urls = [
        "https://github.com/ros/gencpp/archive/0.6.5.tar.gz",
    ],
)

http_archive(
    name = "ros_genpy",
    build_file = "//third_party/ros:genpy.BUILD",
    sha256 = "a0cf129fe90cf342090aac4ca63c33f993f4feecdaf5bc74f410b633e3ea0afc",
    strip_prefix = "genpy-0.6.16",
    urls = [
        "https://github.com/ros/genpy/archive/0.6.16.tar.gz",
    ],
)

http_archive(
    name = "ros_common_msgs",
    build_file = "//third_party/ros:common_msgs.BUILD",
    sha256 = "74af8cc88bdc9c23cbc270d322e50562857e2c877359423f389d51c0735ee230",
    strip_prefix = "common_msgs-1.13.1",
    urls = [
        "https://github.com/ros/common_msgs/archive/1.13.1.tar.gz",
    ],
)

http_archive(
    name = "ros_std_msgs",
    build_file = "//third_party/ros:std_msgs.BUILD",
    sha256 = "ee6592d37b00a94cab8216aac2cfb5120f6da09ffa94bfe197fe8dc76dd21326",
    strip_prefix = "std_msgs-0.5.13",
    urls = [
        "https://github.com/ros/std_msgs/archive/0.5.13.tar.gz",
    ],
)

http_archive(
    name = "ros_comm",
    build_file = "//third_party/ros/ros_comm:ros_comm.BUILD",
    patch_args = ["-p1"],
    patches = [
        "//third_party/ros/ros_comm:0001-fix-include-path.patch",
        "//third_party/ros/ros_comm:0002-remove-pluginlib.patch",
    ],
    sha256 = "0a51857a50cf646db4af85469cb0e4877b1484f7aa0c00ec65a8be7ff574a886",
    strip_prefix = "ros_comm-1.16.0",
    urls = [
        "https://github.com/ros/ros_comm/archive/refs/tags/1.16.0.tar.gz",
    ],
)

http_archive(
    name = "rosconsole",
    build_file = "//third_party/ros:rosconsole.BUILD",
    sha256 = "0b2cbc4f9a92466c0fbae7863482b286ef87692de4941527cb429e6c74639246",
    strip_prefix = "rosconsole-1.14.3",
    urls = [
        "https://github.com/ros/rosconsole/archive/1.14.3.tar.gz",
    ],
)

http_archive(
    name = "ros_comm_msgs",
    build_file = "//third_party/ros:ros_comm_msgs.BUILD",
    sha256 = "5b8b91e8671d03ea84ba32a3ea7360bc4594655e7ba3ec6677a984f393aaafbd",
    strip_prefix = "ros_comm_msgs-1.11.3",
    urls = [
        "https://github.com/ros/ros_comm_msgs/archive/1.11.3.tar.gz",
    ],
)

# setup compile toolchain.
BAZEL_TOOLCHAIN_TAG = "0.8.2"

BAZEL_TOOLCHAIN_SHA = "0fc3a2b0c9c929920f4bed8f2b446a8274cad41f5ee823fd3faa0d7641f20db0"

http_archive(
    name = "com_grail_bazel_toolchain",
    canonical_id = BAZEL_TOOLCHAIN_TAG,
    sha256 = BAZEL_TOOLCHAIN_SHA,
    strip_prefix = "bazel-toolchain-{tag}".format(tag = BAZEL_TOOLCHAIN_TAG),
    url = "https://github.com/grailbio/bazel-toolchain/archive/refs/tags/{tag}.tar.gz".format(tag = BAZEL_TOOLCHAIN_TAG),
)

load("@com_grail_bazel_toolchain//toolchain:deps.bzl", "bazel_toolchain_dependencies")

bazel_toolchain_dependencies()

load("@com_grail_bazel_toolchain//toolchain:rules.bzl", "llvm_toolchain")

llvm_toolchain(
    name = "llvm_toolchain",
    llvm_versions = {
        "": "15.0.6",
        "darwin-aarch64": "15.0.7",
        "darwin-x86_64": "15.0.7",
    },
    sha256 = {
        "": "38bc7f5563642e73e69ac5626724e206d6d539fbef653541b34cae0ba9c3f036",
        "darwin-aarch64": "867c6afd41158c132ef05a8f1ddaecf476a26b91c85def8e124414f9a9ba188d",
        "darwin-x86_64": "d16b6d536364c5bec6583d12dd7e6cf841b9f508c4430d9ee886726bd9983f1c",
    },
    strip_prefix = {
        "": "clang+llvm-15.0.6-x86_64-linux-gnu-ubuntu-18.04",
        "darwin-aarch64": "clang+llvm-15.0.7-arm64-apple-darwin22.0",
        "darwin-x86_64": "clang+llvm-15.0.7-x86_64-apple-darwin21.0",
    },
    urls = {
        "": ["https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.6/clang+llvm-15.0.6-x86_64-linux-gnu-ubuntu-18.04.tar.xz"],
        "darwin-aarch64": ["https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/clang+llvm-15.0.7-arm64-apple-darwin22.0.tar.xz"],
        "darwin-x86_64": ["https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/clang+llvm-15.0.7-x86_64-apple-darwin21.0.tar.xz"],
    },
)

load("@llvm_toolchain//:toolchains.bzl", "llvm_register_toolchains")

llvm_register_toolchains()
