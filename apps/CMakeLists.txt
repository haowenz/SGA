add_executable(sga_example example.cc)
target_link_libraries(sga_example PRIVATE sga)
add_executable(txt2gfa txt2gfa.cc)
target_link_libraries(txt2gfa PRIVATE sga)
