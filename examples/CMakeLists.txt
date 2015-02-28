set(EXAMPLES
    ComputeMatrix
    TavisCummings
    Vectors
    OpenMP
   )

set(EXPENSIVE_EXAMPLES
    MatVecBenchmark
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

if(MBO_BUILD_DOCS)
  set(EXAMPLE MatVecBenchmark)
  set(DOXIFIED_SOURCE ${PROJECT_BINARY_DIR}/doc/${EXAMPLE}.dox.c)
  add_custom_command(OUTPUT ${DOXIFIED_SOURCE}
      COMMAND runhaskell
      ${PROJECT_SOURCE_DIR}/utilities/ProgramToDoxygen.hs
      ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}.c
      ${EXAMPLE}
      mbo_benchmarks ${DOXIFIED_SOURCE}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}.c)
  add_custom_target(${EXAMPLE}-doxygen DEPENDS ${DOXIFIED_SOURCE})
  add_dependencies(doc ${EXAMPLE}-doxygen)
endif()
