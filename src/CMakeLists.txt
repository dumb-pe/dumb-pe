cmake_minimum_required(VERSION 2.8)

project(eos)

set(CMAKE_C_FLAGS "-std=c99 -pthread -I/usr/local/include/guile/2.2")

file(GLOB SOURCES "*.c")
set(GUILE_INCLUDE "/usr/local/include/guile/2.2")
set(GUILE_LIBRARY "/usr/local/lib/libguile-2.2.so")
include_directories(${CMAKE_SOURCE_DIR} ${GUILE_INCLUDE})

set(xfiles)
foreach(_file ${SOURCES})
  string(REPLACE ".c" ".x" file_x ${_file})
  add_custom_command(
    OUTPUT ${file_x}
    COMMAND guile-snarf -o ${file_x} ${_file} -pthread -I/usr/local/include/guile/2.2
    DEPENDS ${_file}
  )
  list(APPEND xfiles ${file_x})
endforeach()
add_custom_target(do_snarf DEPENDS ${file_x})

add_executable(test-eos ${SOURCES})
target_link_libraries(test-eos ${GUILE_LIBRARY} m lapack blas)
add_dependencies(test-eos do_snarf)
