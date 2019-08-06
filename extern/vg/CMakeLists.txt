protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS "vg.proto")
add_library(vg PUBLIC "${PROTO_SRCS}")
