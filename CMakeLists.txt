cmake_minimum_required(VERSION 3.3)

project(SGA DESCRIPTION "Implementation of sequence to graph alignment algorithms")

find_package(Protobuf REQUIRED)
find_package(OpenMP)

add_subdirectory(src)

option(SGA_BUILD_TESTING "Build the tests")
if(SGA_BUILD_TESTING)
  find_package(GTest REQUIRED)
  enable_testing()
  add_subdirectory(tests)
endif()
