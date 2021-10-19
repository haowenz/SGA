cmake_minimum_required(VERSION 3.3)

project(SGA DESCRIPTION "Implementation of sequence to graph alignment algorithms")

add_subdirectory(src)
add_subdirectory(apps)

option(SGA_BUILD_TESTING "Build the tests")
if(SGA_BUILD_TESTING)
  message(STATUS "Build the tests.")
  # GoogleTest requires at least C++11
  set(CMAKE_CXX_STANDARD 11)

  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_subdirectory(tests)
endif()
