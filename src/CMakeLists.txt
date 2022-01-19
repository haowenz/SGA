find_package(ZLIB REQUIRED)

include(ExternalProject)
ExternalProject_Add(gfatools
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/gfatools"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E env
  HOST_CC=${CMAKE_CXX_COMPILER}
  CC=${CMAKE_CXX_COMPILER}
  make -C <SOURCE_DIR> libgfa1.a
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND ""
  )

add_library(sga INTERFACE)
target_include_directories(sga INTERFACE "${CMAKE_SOURCE_DIR}/include/sga" "${CMAKE_SOURCE_DIR}/extern/klib" "${CMAKE_SOURCE_DIR}/extern/gfatools")
target_link_libraries(sga INTERFACE ZLIB::ZLIB "${CMAKE_SOURCE_DIR}/extern/gfatools/libgfa1.a")
add_dependencies(sga gfatools)
