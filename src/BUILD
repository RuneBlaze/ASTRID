package(default_visibility = ["//test:__subpackages__"])
load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
load("@pybind11_bazel//:build_defs.bzl", "pybind_library")


cc_binary(
    name = "ASTRID-phydstar",
    srcs = [
        "astrid.cpp",
        "@bazel_tools//tools/jdk",
    ],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods",
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_binary(
    name = "ASTRID",
    srcs = [
        "astrid.cpp",
    ],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods:nophydstar",
        "@pybind11//:pybind11_embed",
        "@phylokit",
    ],
)

cc_library(
    name = "Args",
    srcs = ["help.hpp"],
    hdrs = ["Args.hpp"],
)

pybind_extension(
    name="asterid",
    srcs = ["astrid.cpp"],
    deps = ["astrid_base"],
)

cc_library(
    name = "astrid_base",
    # srcs = ["astrid.cpp"],
    hdrs = ["astrid.hpp"],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods:nophydstar",
        "@pybind11//:pybind11_embed",
        "@phylokit",
    ],
)

# pybind_library(
#     name = "astrid-python",
#     srcs = ["pythonwrapper.cpp"],
#     deps = [
#         "astrid_base",
#     ],
# )

cc_library(
    name = "Multind",
    srcs = ["multind.cpp"],
    hdrs = ["multind.hpp"],
    deps = [
        "@boost//:tokenizer",
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_library(
    name = "Octal",
    srcs = ["octal.cpp"],
    hdrs = ["octal.hpp"],
    deps = [
        "@phylokit",
    ],
)

