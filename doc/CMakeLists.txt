set(EXAMPLES
    SimpleExample
   )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    )

foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.c)
  target_link_libraries(${e} MBO)
  if(RUN_EXAMPLES_AS_TESTS)
    add_test(${e} ${e})
  endif()
endforeach()

