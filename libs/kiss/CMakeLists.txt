cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(kiss STATIC
    src/logger
    src/path
    src/string
    src/uuid
)

SET_TARGET_PROPERTIES( kiss PROPERTIES COMPILE_FLAGS -fPIC)
