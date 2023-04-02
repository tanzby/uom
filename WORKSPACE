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
    build_file = "//third_party:eigen.BUILD",
    remote = "https://gitlab.com/libeigen/eigen.git",
    tag = "3.4.0",
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
    patch_cmds = [
        "mv support/bazel/.bazelrc .bazelrc",
        "mv support/bazel/.bazelversion .bazelversion",
        "mv support/bazel/BUILD.bazel BUILD.bazel",
        "mv support/bazel/WORKSPACE.bazel WORKSPACE.bazel",
    ],
    remote = "https://github.com/fmtlib/fmt.git",
    tag = "9.1.0",
)

git_repository(
    name = "bazel_skylib",
    remote = "https://github.com/bazelbuild/bazel-skylib.git",
    tag = "1.4.1",
)

git_repository(
    name = "gflags",
    remote = "https://github.com/gflags/gflags.git",
    tag = "v2.2.2",
)

git_repository(
    name = "glog",
    patch_args = ["-p1"],
    patches = [
        "//third_party:glog.patch",
    ],
    remote = "https://github.com/google/glog.git",
    repo_mapping = {
        "@com_github_gflags_gflags": "@gflags",
    },
    tag = "v0.6.0",
)

git_repository(
    name = "com_google_protobuf",
    remote = "https://github.com/google/protobuf.git",
    tag = "v3.19.6",
)

git_repository(
    name = "rules_python",
    remote = "https://github.com/bazelbuild/rules_python.git",
    tag = "0.20.0",
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
    name = "boost",
    build_file = "//third_party:boost.BUILD",
    patch_args = ["-p1"],
    patches = ["//third_party:boost.patch"],
    sha256 = "273f1be93238a068aba4f9735a4a2b003019af067b9c183ed227780b8f36062c",
    strip_prefix = "boost_1_79_0",
    urls = [
        "https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz",
    ],
)

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
    sha256 = "a364f5162c7d1a455cc915e8e3cf5f4bd8b75d09bc0f53965b0c9ca1383c52c8",
    strip_prefix = "zstd-1.4.4",
    urls = [
        "https://github.com/facebook/zstd/archive/v1.4.4.tar.gz",
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

############################################ ROS ############################################
# all packages' version follow the *melodic* release. see https://index.ros.org/

http_archive(
    name = "com_github_mvukov_rules_ros",
    patch_args = ["-p1"],
    patches = ["//third_party/ros:com_github_mvukov_rules_ros.patch"],
    sha256 = "d0c0424df79d5a9a7f375cdcaf542997c72937e276a494019d58b95038b1fec5",
    strip_prefix = "rules_ros-0.1.0",
    urls = [
        "https://github.com/mvukov/rules_ros/archive/refs/tags/v0.1.0.zip",
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
    sha256 = "d5a0ad09fa878d9f3d6d7f3e8c7854f0f160aeeea9c4d332e3dc87552087ca68",
    # https://index.ros.org/p/roscpp_core/github-ros-roscpp_core/#melodic
    strip_prefix = "roscpp_core-0.6.14",
    urls = [
        "https://github.com/ros/roscpp_core/archive/0.6.14.tar.gz",
    ],
)

http_archive(
    name = "ros_genmsg",
    build_file = "//third_party/ros:genmsg.BUILD",
    sha256 = "0e414846823a2aaa7781f81268251c7c9a45ff96cef8e6a78bbbbcf7e4c28d56",
    # https://index.ros.org/p/genmsg/github-ros-genmsg/#melodic
    strip_prefix = "genmsg-0.5.16",
    urls = [
        "https://github.com/ros/genmsg/archive/0.5.16.tar.gz",
    ],
)

http_archive(
    name = "ros_gencpp",
    build_file = "//third_party/ros:gencpp.BUILD",
    sha256 = "05acfeeb1bbc374356bf7674fee2a7aab3bf6a48ebad4a06fd0f0d4455a60720",
    # https://index.ros.org/p/gencpp/github-ros-gencpp/#melodic
    strip_prefix = "gencpp-0.6.5",
    urls = [
        "https://github.com/ros/gencpp/archive/0.6.5.tar.gz",
    ],
)

http_archive(
    name = "ros_genpy",
    build_file = "//third_party/ros:genpy.BUILD",
    sha256 = "a0cf129fe90cf342090aac4ca63c33f993f4feecdaf5bc74f410b633e3ea0afc",
    # https://index.ros.org/p/genpy/github-ros-genpy/#melodic
    strip_prefix = "genpy-0.6.16",
    urls = [
        "https://github.com/ros/genpy/archive/0.6.16.tar.gz",
    ],
)

http_archive(
    name = "ros_common_msgs",
    build_file = "//third_party/ros:common_msgs.BUILD",
    sha256 = "d3c698e7c164e97d135ca341420f266b95f5c4ef995e8faf9a5362f44987718d",
    # https://index.ros.org/r/common_msgs/github-ros-common_msgs/#melodic
    strip_prefix = "common_msgs-1.12.8",
    urls = [
        "https://github.com/ros/common_msgs/archive/1.12.8.tar.gz",
    ],
)

http_archive(
    name = "ros_std_msgs",
    build_file = "//third_party/ros:std_msgs.BUILD",
    sha256 = "ee6592d37b00a94cab8216aac2cfb5120f6da09ffa94bfe197fe8dc76dd21326",
    # https://index.ros.org/p/std_msgs/github-ros-std_msgs/#melodic
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
        "//third_party/ros/ros_comm:0002-fix-boost-placeholders-namespace.patch",
        "//third_party/ros/ros_comm:0003-remove-pluginlib.patch",
    ],
    sha256 = "1ab904142eecc0be8fd09662df7411d4494bff765dd2a2c6508174802ca7e215",
    # https://index.ros.org/p/ros_comm/github-ros-ros_comm/#melodic
    strip_prefix = "ros_comm-1.14.13",
    urls = [
        "https://github.com/ros/ros_comm/archive/refs/tags/1.14.13.tar.gz",
    ],
)

http_archive(
    name = "rosconsole",
    build_file = "//third_party/ros:rosconsole.BUILD",
    sha256 = "234d83dfddcf864e5d223eaedd58e1505ad0d2707ea4ff497b69c4f28501f179",
    # https://index.ros.org/p/rosconsole/github-ros-rosconsole/#melodic
    strip_prefix = "rosconsole-1.13.18",
    urls = [
        "https://github.com/ros/rosconsole/archive/1.13.18.tar.gz",
    ],
)

http_archive(
    name = "ros_comm_msgs",
    build_file = "//third_party/ros:ros_comm_msgs.BUILD",
    sha256 = "5b8b91e8671d03ea84ba32a3ea7360bc4594655e7ba3ec6677a984f393aaafbd",
    # https://index.ros.org/r/ros_comm_msgs/github-ros-ros_comm_msgs/#melodic
    strip_prefix = "ros_comm_msgs-1.11.3",
    urls = [
        "https://github.com/ros/ros_comm_msgs/archive/1.11.3.tar.gz",
    ],
)
