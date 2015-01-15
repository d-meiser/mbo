set(EXAMPLES
    ComputeMatrix
    TavisCummings
    Vectors
    OpenMP
   )

set(EXPENSIVE_EXAMPLES
    Benchmark
   )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    )

foreach(e ${EXAMPLES} ${EXPENSIVE_EXAMPLES})
  add_executable(${e} ${e}.c)
  target_link_libraries(${e} MBO)
endforeach()

if(RUN_EXAMPLES_AS_TESTS)
  foreach(e ${EXAMPLES})
    add_test(${e} ${e})
  endforeach()
endif()
