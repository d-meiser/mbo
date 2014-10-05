include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
  ${PROJECT_BINARY_DIR}
  ${GTEST_INCLUDE_DIRS}
  )

set(TESTS
  test_DenseMatrices
  test_Embedding
  test_MboElemOp
  test_MboNumOp
  test_MboNumSubMatrix
  test_MboProdSpace
  test_MboTensorOp
  test_MboVec
  test_SimpleTOps
  test_Tile
  test_SmallOps
  )

foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} MBO gtest_main gtest)
  if(MBO_ENABLE_COVERAGE)
    target_link_libraries(${t} "${COV_LIBRARIES}")
  endif()
  add_test(${t} ${t})
endforeach()

