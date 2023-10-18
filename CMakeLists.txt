cmake_minimum_required(VERSION 2.4)
project(splines-lib)

file(GLOB SOURCES
    "*.h"
    "*.cpp")
include_directories(.)
add_library(splines-lib STATIC ${SOURCES})
