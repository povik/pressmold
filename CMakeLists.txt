cmake_minimum_required (VERSION 3.9)
project(pressmold)

set(OPENSTA_HOME ${PROJECT_SOURCE_DIR}/third_party/OpensTA)
add_subdirectory(third_party/OpenSTA)
include(src/cmake/swig_lib.cmake)
add_subdirectory(src)
