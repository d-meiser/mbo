cmake_minimum_required(VERSION 2.8)
project(QUO)

enable_testing()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=auto -Wall")
  set(CMAKE_C_FLAGS_DEBUG 
      "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG
      "${CMAKE_EXE_LINKER_FLAGS} -lgcov")
endif()

add_subdirectory(include)
add_subdirectory(src)
