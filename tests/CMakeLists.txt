include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

set(TESTS
  test_MboElemOp
  test_MboProdSpace
  test_MboTensorOp
  test_MboVec
  )

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} MBO)
  add_test(${t} ${t})
endforeach()

