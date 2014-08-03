include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
set(QUO_SRCS
  MboElemOp.c
  MboProdSpace.c
  MboTensorOp.c
  )

add_library(QUO
    ${QUO_SRCS}
    )
if(UNIX)
  target_link_libraries(QUO m)
endif()

