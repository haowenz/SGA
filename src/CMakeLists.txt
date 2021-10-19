add_library(sga INTERFACE)
target_include_directories(sga INTERFACE "${CMAKE_SOURCE_DIR}/include/sga" "${CMAKE_SOURCE_DIR}/extern/klib")
