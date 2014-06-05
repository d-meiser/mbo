cmake_minimum_required(VERSION 2.8)
project(QUO)

enable_testing()

option(RUN_EXAMPLES_AS_TESTS 
    "Whether to run the examples as part of the test suite" ON)

if(CMAKE_COMPILER_IS_GNUCC)
  set(WARNING_OPTIONS "-Wall -Wpedantic")
  set(CMAKE_C_FLAGS 
      "${CMAKE_C_FLAGS} -std=c99 -fdiagnostics-color=auto ${WARNING_OPTIONS}")
  set(CMAKE_C_FLAGS_DEBUG 
      "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG
      "${CMAKE_EXE_LINKER_FLAGS} -lgcov")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  set(WARNING_OPTIONS "-Wall -Wpedantic")
  set(CMAKE_C_FLAGS 
      "${CMAKE_C_FLAGS} -std=c99 -fdiagnostics-color=auto ${WARNING_OPTIONS}")
  set(CMAKE_C_FLAGS_DEBUG 
      "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG
      "${CMAKE_EXE_LINKER_FLAGS} -lgcov")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
endif()

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(examples)
