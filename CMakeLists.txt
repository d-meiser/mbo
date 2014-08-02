cmake_minimum_required(VERSION 2.8)
project(QUO)
enable_testing()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
include(AddCCompilerFlagIfSupported)

option(RUN_EXAMPLES_AS_TESTS 
    "Whether to run the examples as part of the test suite" ON)

if(CMAKE_COMPILER_IS_GNUCC)
  set(WARNING_FLAGS "-Wall")
  add_c_compiler_flag_if_supported("-Wpedantic")
  add_c_compiler_flag_if_supported("-fdiagnostics-color=auto")
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_c_compiler_flag_if_supported("-fsanitize=address")
    add_c_compiler_flag_if_supported("-ftest-coverage -fprofile-arcs")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${CC_FLAGS} ${WARNING_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  set(WARNING_FLAGS "-Wall -Wpedantic")
  set(CMAKE_C_FLAGS 
      "${CMAKE_C_FLAGS} -std=c99 -fdiagnostics-color=auto ${WARNING_FLAGS}")
  set(CMAKE_C_FLAGS_DEBUG 
      "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -ftest-coverage -fprofile-arcs")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG
      "${CMAKE_EXE_LINKER_FLAGS} -ftest-coverage -fprofile-arcs")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
endif()

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(examples)
