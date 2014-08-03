include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

set(TESTS
  test_ElemOp
  test_ProdSpace
  test_TensorOp
  )

if(CMAKE_COMPILER_IS_GNUCC)
  if(MBO_ENABLE_COVERAGE)
    add_c_compiler_flag_if_supported("-ftest-coverage -fprofile-arcs")
    add_c_linker_flag_if_supported("-ftest-coverage -fprofile-arcs")
  endif()
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  if(MBO_ENABLE_COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftest-coverage -fprofile-arcs")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -ftest-coverage -fprofile-arcs")
  endif()
endif()

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} QUO)
  add_test(${t} ${t})
endforeach()

