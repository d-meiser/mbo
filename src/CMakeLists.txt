include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
set(QUO_SRCS
  ElemOp.c
  ProdSpace.c
  TensorOp.c
  )

add_library(QUO
    ${QUO_SRCS}
    )
if(UNIX)
  target_link_libraries(QUO m)
endif()

set(TESTS
  test_ElemOp
  test_ProdSpace
  test_TensorOp
  )

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} QUO)
  add_test(${t} ${t})
endforeach()
