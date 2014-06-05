set(EXAMPLES
    TavisCummings
   )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    )

foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.c)
  target_link_libraries(${e} QUO)
  if(RUN_EXAMPLES_AS_TESTS)
  add_test(${e} ${e})
  endif()
endforeach()

