include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${PROJECT_BINARY_DIR}
  ${GTEST_INCLUDE_DIRS}
  )

set(TESTS
  test_MboElemOp
  test_MboProdSpace
  test_MboTensorOp
  test_MboVec
  )

foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} MBO gtest_main gtest)
  if(MBO_ENABLE_COVERAGE)
    target_link_libraries(${t} "${COV_LIBRARIES}")
  endif()
  add_test(${t} ${t})
endforeach()

