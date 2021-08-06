cc_library(
    name = "DistanceMethods",
    hdrs = ["DistanceMethods.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "PhyDStar",
        ":FastMEInterface",
        
        ":UPGMA",
    ],
)

cc_library(
    name = "nophydstar",
    hdrs = ["DistanceMethods.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "PhyDStar-stub",
        ":FastMEInterface",
        
        ":UPGMA",
    ],
)

cc_library(
    name = "FastMEInterface",
    srcs = ["FastMEInterface.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        "//third_party/fastme",
        "@phylokit",
    ],
)

cc_library(
    name = "RapidNJInterface",
    srcs = ["RapidNJInterface.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps =
        select(
            {
                
                "//conditions:default": [
                    
                ],
            },
        ) + [
            "@phylokit",
        ],
)

cc_library(
    name = "PhyDStar-stub",
    srcs = [
        "phydstar.cpp",
    ],
    hdrs = ["DistanceMethods.hpp"],
    deps = ["@phylokit"],
)

cc_library(
    name = "PhyDStar",
    srcs = [
        "phydstar.cpp",
    ],
    hdrs = ["DistanceMethods.hpp"],
    copts = ["-DUSE_PHYDSTAR"],
    deps = [
        ":whereami",
        "//:jni_headers",
        "@phylokit",
    ],
    linkopts =["-ldl"]
)

cc_library(
    name = "UPGMA",
    srcs = ["UPGMA.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_library(
    name = "whereami",
    srcs = ["whereami.c"],
    hdrs = ["whereami.h"],
    copts=["-std=c99"]
)
