licenses(["notice"])  # Boost Software 1.0

load("@uom//third_party:boost.bzl", "boost_library")

package(default_visibility = ["//visibility:public"])

boost_library(
    name = "algorithm",
    deps = [
        ":array",
        ":assert",
        ":bind",
        ":concept_check",
        ":config",
        ":detail",
        ":exception",
        ":function",
        ":integer",
        ":iterator",
        ":mpl",
        ":range",
        ":ref",
        ":regex",
        ":static_assert",
        ":string",
        ":throw_exception",
        ":tuple",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "bimap",
)

boost_library(
    name = "concept_archetype",
)

boost_library(
    name = "implicit_cast",
)

boost_library(
    name = "indirect_reference",
)

boost_library(
    name = "iterator_adaptors",
)

boost_library(
    name = "lambda",
)

boost_library(
    name = "make_shared",
)

boost_library(
    name = "phoenix",
)

boost_library(
    name = "pointer_to_other",
)

boost_library(
    name = "pool",
)

boost_library(
    name = "property_map",
)

boost_library(
    name = "proto",
)

boost_library(
    name = "shared_array",
)

boost_library(
    name = "weak_ptr",
)

boost_library(
    name = "xpressive",
)

boost_library(
    name = "graph",
    deps = [
        ":algorithm",
        ":any",
        ":assert",
        ":concept_archetype",
        ":core",
        ":detail",
        ":enable_shared_from_this",
        ":foreach",
        ":function",
        ":implicit_cast",
        ":indirect_reference",
        ":intrusive_ptr",
        ":iterator_adaptors",
        ":lexical_cast",
        ":make_shared",
        ":mpl",
        ":multi_index",
        ":parameter",
        ":pending",
        ":pointer_to_other",
        ":property_map",
        ":property_tree",
        ":proto",
        ":range",
        ":ref",
        ":scoped_array",
        ":shared_array",
        ":throw_exception",
        ":tuple",
        ":unordered",
        ":utility",
        ":weak_ptr",
        ":xpressive",
    ],
)

boost_library(
    name = "crc",
    deps = [
        ":array",
        ":config",
        ":cstdint",
        ":integer",
    ],
)

boost_library(
    name = "string",
)

boost_library(
    name = "align",
    deps = [
        ":assert",
        ":config",
        ":core",
        ":cstdint",
        ":static_assert",
        ":throw_exception",
    ],
)

boost_library(
    name = "assert",
    deps = [
        ":cstdint",
        ":current_function",
    ],
)

boost_library(
    name = "memory_order",
)

boost_library(
    name = "atomic",
    exclude_srcs = [
        "libs/atomic/src/wait_on_address.cpp",
        "libs/atomic/src/find_address_sse2.cpp",
        "libs/atomic/src/find_address_sse41.cpp",
    ],
    deps = [
        ":align",
        ":assert",
        ":config",
        ":cstdint",
        ":memory_order",
        ":predef",
        ":preprocessor",
        ":static_assert",
        ":type_traits",
    ],
)

boost_library(
    name = "get_pointer",
)

boost_library(
    name = "is_placeholder",
)

boost_library(
    name = "visit_each",
)

boost_library(
    name = "bind",
    deps = [
        ":config",
        ":core",
        ":detail",
        ":get_pointer",
        ":is_placeholder",
        ":mem_fn",
        ":ref",
        ":static_assert",
        ":type",
        ":visit_each",
    ],
)

boost_library(
    name = "version",
)

boost_library(
    name = "operators",
    deps = [
        ":iterator",
    ],
)

boost_library(
    name = "chrono",
    deps = [
        ":config",
        ":mpl",
        ":operators",
        ":predef",
        ":ratio",
        ":system",
        ":throw_exception",
        ":version",
    ],
)

boost_library(
    name = "concept",
)

boost_library(
    name = "concept_check",
    deps = [
        ":concept",
    ],
)

boost_library(
    name = "config",
)

boost_library(
    name = "container",
    deps = [
        ":core",
        ":intrusive",
        ":move",
    ],
)

boost_library(
    name = "circular_buffer",
    deps = [
        ":call_traits",
        ":concept_check",
        ":container",
        ":iterator",
        ":throw_exception",
    ],
)

boost_library(
    name = "intrusive_ptr",
    deps = [
        ":smart_ptr",
    ],
)

boost_library(
    name = "conversion",
)

boost_library(
    name = "core",
    deps = [
        ":config",
    ],
)

boost_library(
    name = "token_functions",
    deps = [
        ":assert",
        ":config",
        ":detail",
        ":mpl",
        ":throw_exception",
    ],
)

boost_library(
    name = "token_iterator",
    deps = [
        ":assert",
        ":iterator",
        ":token_functions",
    ],
)

boost_library(
    name = "tokenizer",
    deps = [
        ":token_iterator",
    ],
)

boost_library(
    name = "date_time",
    deps = [
        ":algorithm",
        ":assert",
        ":config",
        ":detail",
        ":io",
        ":lexical_cast",
        ":limits",
        ":mpl",
        ":operators",
        ":range",
        ":serialization",
        ":shared_ptr",
        ":static_assert",
        ":throw_exception",
        ":tokenizer",
        ":type_traits",
    ],
)

boost_library(
    name = "detail",
    deps = [
        ":limits",
    ],
)

boost_library(
    name = "exception",
    deps = [
        ":config",
    ],
)

boost_library(
    name = "filesystem",
    copts = [
        "-DBOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF",
        "-Wno-unused-function",
        "-Wno-unused-parameter",
        "-Wno-deprecated-declarations",
    ],
    deps = [
        ":atomic",
        ":config",
        ":functional",
        ":io",
        ":iterator",
        ":range",
        ":scoped_array",
        ":shared_ptr",
        ":smart_ptr",
        ":system",
        ":type_traits",
    ],
)

boost_library(
    name = "foreach",
)

boost_library(
    name = "function_equal",
)

boost_library(
    name = "mem_fn",
)

boost_library(
    name = "typeof",
    deps = [
        ":config",
        ":mpl",
        ":preprocessor",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "function",
    deps = [
        ":assert",
        ":bind",
        ":config",
        ":detail",
        ":function_equal",
        ":integer",
        ":mem_fn",
        ":mpl",
        ":ref",
        ":type_index",
        ":typeof",
    ],
)

boost_library(
    name = "function_types",
    deps = [
        ":blank",
    ],
)

boost_library(
    name = "rational",
    deps = [
        ":operators",
    ],
)

boost_library(
    name = "functional",
    deps = [
        ":assert",
        ":container_hash",
        ":cstdint",
        ":detail",
        ":integer",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "geometry",
    deps = [
        ":algorithm",
        ":call_traits",
        ":function_types",
        ":lexical_cast",
        ":multiprecision",
        ":numeric",
        ":qvm",
        ":range",
        ":rational",
        ":tokenizer",
        ":tuple",
        ":variant",
        ":version",
    ],
)

boost_library(
    name = "container_hash",
)

boost_library(
    name = "heap",
    deps = [
        ":parameter",
    ],
)

boost_library(
    name = "integer",
    deps = [
        ":integer_traits",
        ":static_assert",
    ],
)

boost_library(
    name = "iterator",
    deps = [
        ":detail",
        ":mpl",
        ":static_assert",
        ":type_traits",
    ],
)

boost_library(
    name = "intrusive",
    deps = [
        ":assert",
        ":static_assert",
    ],
)

boost_library(
    name = "io",
)

boost_library(
    name = "ref",
)

boost_library(
    name = "cstdint",
)

boost_library(
    name = "iostreams",
    copts = [
        "-Wno-constant-conversion",
    ],
    exclude_srcs = [
        "libs/iostreams/src/bzip2.cpp",
        "libs/iostreams/src/lzma.cpp",
        "libs/iostreams/src/zstd.cpp",
    ],
    deps = [
        ":assert",
        ":bind",
        ":call_traits",
        ":config",
        ":cstdint",
        ":detail",
        ":integer_traits",
        ":numeric",
        ":range",
        ":ref",
        ":shared_ptr",
        ":static_assert",
        ":throw_exception",
        ":type_traits",
        "@zlib",
    ],
)

boost_library(
    name = "pending",
    deps = [
        ":config",
        ":integer",
        ":iterator",
        ":limits",
        ":mpl",
        ":next_prior",
        ":none",
        ":operators",
        ":optional",
        ":serialization",
        ":static_assert",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "random",
    deps = [
        ":assert",
        ":config",
        ":cstdint",
        ":detail",
        ":integer",
        ":integer_traits",
        ":limits",
        ":mpl",
        ":noncopyable",
        ":operators",
        ":pending",
        ":range",
        ":static_assert",
        ":system",
        ":throw_exception",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "format",
    deps = [
        ":assert",
        ":config",
        ":detail",
        ":limits",
        ":optional",
        ":shared_ptr",
        ":throw_exception",
        ":utility",
    ],
)

boost_library(
    name = "math",
    deps = [
        ":format",
    ],
)

boost_library(
    name = "move",
    deps = [
        ":assert",
        ":config",
        ":core",
        ":static_assert",
    ],
)

boost_library(
    name = "mpl",
    deps = [
        ":detail",
        ":move",
        ":preprocessor",
    ],
)

boost_library(
    name = "multi_index",
    hdrs = [
        "boost/multi_index_container.hpp",
        "boost/multi_index_container_fwd.hpp",
    ],
    deps = [
        ":foreach",
        ":serialization",
        ":static_assert",
        ":tuple",
    ],
)

boost_library(
    name = "type",
)

boost_library(
    name = "none",
    deps = [
        ":config",
        ":none_t",
    ],
)

boost_library(
    name = "none_t",
)

boost_library(
    name = "optional",
    deps = [
        ":assert",
        ":config",
        ":core",
        ":detail",
        ":move",
        ":mpl",
        ":none",
        ":static_assert",
        ":throw_exception",
        ":type",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "parameter",
    deps = [
        ":mp11",
    ],
)

boost_library(
    name = "predef",
)

boost_library(
    name = "preprocessor",
)

boost_library(
    name = "range",
    deps = [
        ":concept_check",
        ":cstdint",
        ":iterator",
        ":optional",
        ":tuple",
    ],
)

boost_library(
    name = "integer_traits",
)

boost_library(
    name = "ratio",
    deps = [
        ":cstdint",
        ":integer_traits",
        ":type_traits",
    ],
)

boost_library(
    name = "limits",
)

boost_library(
    name = "scoped_array",
)

boost_library(
    name = "cregex",
)

boost_library(
    name = "regex",
    copts = [
        "-Wno-deprecated-register",
        "-Wno-register",
    ],
    deps = [
        ":assert",
        ":config",
        ":cregex",
        ":cstdint",
        ":functional",
        ":integer",
        ":limits",
        ":mpl",
        ":ref",
        ":scoped_array",
        ":scoped_ptr",
        ":shared_ptr",
        ":smart_ptr",
        ":throw_exception",
        ":type_traits",
    ],
)

boost_library(
    name = "compressed_pair",
)

boost_library(
    name = "pointee",
)

boost_library(
    name = "spirit",
    deps = [
        ":compressed_pair",
        ":endian",
        ":optional",
        ":ref",
    ],
)

boost_library(
    name = "archive",
    deps = [
        ":assert",
        ":cstdint",
        ":integer",
        ":integer_traits",
        ":io",
        ":noncopyable",
        ":pointee",
        ":preprocessor",
        ":scoped_ptr",
        ":spirit",
        ":static_assert",
    ],
)

boost_library(
    name = "call_traits",
)

boost_library(
    name = "swap",
)

boost_library(
    name = "array",
    deps = [
        ":functional",
        ":swap",
        ":throw_exception",
    ],
)

boost_library(
    name = "shared_ptr",
    deps = [
        ":assert",
        ":smart_ptr",
    ],
)

boost_library(
    name = "serialization",
    exclude_srcs = [
        "libs/serialization/src/shared_ptr_helper.cpp",
    ],
    deps = [
        ":archive",
        ":array",
        ":call_traits",
        ":config",
        ":function",
        ":integer_traits",
        ":mpl",
        ":operators",
        ":shared_ptr",
        ":type_traits",
    ],
)

boost_library(
    name = "scoped_ptr",
    deps = [
        ":smart_ptr",
    ],
)

boost_library(
    name = "checked_delete",
)

boost_library(
    name = "smart_ptr",
    srcs = [
        "boost/enable_shared_from_this.hpp",
        "boost/intrusive_ptr.hpp",
        "boost/make_shared.hpp",
        "boost/scoped_array.hpp",
        "boost/scoped_ptr.hpp",
        "boost/shared_array.hpp",
        "boost/shared_ptr.hpp",
        "boost/weak_ptr.hpp",
    ],
    deps = [
        ":align",
        ":checked_delete",
        ":core",
        ":cstdint",
        ":predef",
        ":throw_exception",
        ":utility",
    ],
)

boost_library(
    name = "static_assert",
    deps = [
        ":config",
        ":detail",
    ],
)

boost_library(
    name = "noncopyable",
    deps = [
        ":core",
    ],
)

boost_library(
    name = "cerrno",
)

boost_library(
    name = "system",
    deps = [
        ":assert",
        ":cerrno",
        ":config",
        ":core",
        ":cstdint",
        ":noncopyable",
        ":predef",
        ":utility",
        ":variant2",
    ],
)

boost_library(
    name = "enable_shared_from_this",
    deps = [
        ":smart_ptr",
    ],
)

boost_library(
    name = "exception_ptr",
    deps = [
        ":exception",
    ],
)

boost_library(
    name = "thread",
    srcs = [
        "libs/thread/src/pthread/once_atomic.cpp",
        "libs/thread/src/pthread/thread.cpp",
    ],
    deps = [
        ":atomic",
        ":chrono",
        ":core",
        ":date_time",
        ":detail",
        ":enable_shared_from_this",
        ":exception",
        ":exception_ptr",
        ":io",
        ":system",
    ],
    linkopts = ["-lpthread"],
)

boost_library(
    name = "current_function",
)

boost_library(
    name = "throw_exception",
    deps = [
        ":current_function",
        ":exception",
    ],
)

boost_library(
    name = "type_index",
    deps = [
        ":core",
    ],
)

boost_library(
    name = "aligned_storage",
)

boost_library(
    name = "type_traits",
    deps = [
        ":aligned_storage",
        ":core",
        ":static_assert",
        ":utility",
        ":version",
    ],
)

boost_library(
    name = "tuple",
)

boost_library(
    name = "tr1",
    deps = [
        ":config",
    ],
)

boost_library(
    name = "next_prior",
)

boost_library(
    name = "utility",
    deps = [
        ":next_prior",
        ":noncopyable",
    ],
)

boost_library(
    name = "blank",
)

boost_library(
    name = "variant",
    deps = [
        ":blank",
        ":call_traits",
        ":functional",
        ":math",
        ":type_index",
    ],
)

boost_library(
    name = "variant2",
    deps = [
        "mp11",
    ],
)

boost_library(
    name = "asio",
    deps = [
        ":bind",
        ":config",
        ":date_time",
        ":filesystem",
        ":mpl",
        ":regex",
        ":smart_ptr",
        ":static_assert",
        ":throw_exception",
        ":type_traits",
        ":version",
        "@openssl//:crypto",
        "@openssl//:ssl",
    ],
)

boost_library(
    name = "numeric",
    deps = [
        ":fusion",
        ":math",
        ":multi_array",
        ":range",
        ":serialization",
        ":units",
    ],
)

boost_library(
    name = "multi_array",
)

boost_library(
    name = "fusion",
    deps = [
        ":function_types",
    ],
)

boost_library(
    name = "units",
    deps = [
        ":version",
    ],
)

boost_library(
    name = "lexical_cast",
    deps = [
        ":array",
        ":container",
        ":integer",
        ":math",
        ":numeric",
        ":range",
    ],
)

boost_library(
    name = "uuid",
    deps = [
        ":assert",
        ":config",
        ":cstdint",
        ":iterator",
        ":random",
        ":throw_exception",
        ":tti",
        ":type_traits",
    ],
)

boost_library(
    name = "any",
    deps = [
        ":config",
        ":mpl",
        ":static_assert",
        ":throw_exception",
        ":type_index",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "property_tree",
    deps = [
        ":any",
        ":bind",
        ":format",
        ":multi_index",
        ":optional",
        ":range",
        ":ref",
        ":throw_exception",
        ":utility",
    ],
)

boost_library(
    name = "timer",
    deps = [
        ":chrono",
        ":config",
        ":io",
        ":system",
    ],
)

boost_library(
    name = "unordered",
    hdrs = [
        "boost/unordered_map.hpp",
        "boost/unordered_set.hpp",
    ],
    deps = [
        ":assert",
        ":config",
        ":detail",
        ":functional",
        ":iterator",
        ":limits",
        ":move",
        ":preprocessor",
        ":swap",
        ":throw_exception",
        ":tuple",
        ":type_traits",
        ":utility",
    ],
)

boost_library(
    name = "ptr_container",
    deps = [
        ":array",
        ":assert",
        ":checked_delete",
        ":circular_buffer",
        ":compressed_pair",
        ":config",
        ":iterator",
        ":mpl",
        ":pointee",
        ":range",
        ":scoped_array",
        ":serialization",
        ":static_assert",
        ":type_traits",
        ":unordered",
        ":utility",
    ],
)

boost_library(
    name = "dynamic_bitset",
    deps = [
        ":config",
        ":detail",
        ":limits",
        ":move",
        ":pending",
        ":static_assert",
        ":throw_exception",
        ":utility",
    ],
)

boost_library(
    name = "assign",
    deps = [
        ":array",
        ":config",
        ":detail",
        ":mpl",
        ":preprocessor",
        ":ptr_container",
        ":range",
        ":static_assert",
        ":tuple",
        ":type_traits",
    ],
)

boost_library(
    name = "endian",
    deps = [
        ":assert",
        ":config",
        ":conversion",
        ":core",
        ":cstdint",
        ":current_function",
        ":predef",
        ":static_assert",
        ":system",
        ":type_traits",
    ],
)

boost_library(
    name = "locale",
    srcs = glob(
        [
            "libs/locale/src/encoding/*.cpp",
            "libs/locale/src/encoding/*.ipp",
            "libs/locale/src/shared/*.cpp",
            "libs/locale/src/std/*.cpp",
            "libs/locale/src/util/*.cpp",
            "libs/locale/src/posix/*.cpp",
        ],
    ),
    hdrs = glob(
        [
            "libs/locale/src/encoding/*.hpp",
            "libs/locale/src/shared/*.hpp",
            "libs/locale/src/std/*.hpp",
            "libs/locale/src/util/*.hpp",
            "libs/locale/src/util/iconv.hpp",
            "libs/locale/src/posix/*.hpp",
        ],
    ),
    copts = [
        "-DBOOST_LOCALE_WITH_ICONV",
        "-DBOOST_LOCALE_NO_WINAPI_BACKEND",
        "-Wno-deprecated-declarations",
        "-Wno-unused-const-variable",
        "-Wno-unused-function",
    ],
    deps = [
        ":assert",
        ":config",
        ":cstdint",
        ":smart_ptr",
        ":thread",
        ":unordered",
    ],
)

boost_library(
    name = "tti",
    deps = [
        ":function_types",
    ],
)

boost_library(
    name = "mp11",
)

boost_library(
    name = "multiprecision",
)

boost_library(
    name = "qvm",
)

boost_library(
    name = "scope_exit",
)

boost_library(
    name = "function_output_iterator",
)

boost_library(
    name = "signals2",
    deps = [
        ":bind",
        ":function_output_iterator",
        ":iterator",
        ":iterator_adaptors",
        ":parameter",
        ":variant",
    ],
)

boost_library(
    name = "program_options",
)
