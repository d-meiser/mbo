include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
set(QUO_SRCS
  ElemOp.c
  ProdSpace.c
  TensProdOp.c
  MboList.c
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
  test_TensProdOp
  test_MboList
  )

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} QUO)
  add_test(${t} ${t})
endforeach()
