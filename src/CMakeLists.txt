include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    )

set(MBO_SRCS
  MboElemOp.c
  MboNumOp.c
  MboProdSpace.c
  MboTensorOp.c
  MboTensorOpPrivate.c
  MboVec.c
  Embedding.c
  SimpleTOp.c
  Utilities.c
  )

add_compiler_export_flags()
add_library(MBO
    ${MBO_SRCS}
    )
if(UNIX)
  target_link_libraries(MBO m)
  if(MBO_ENABLE_COVERAGE)
    target_link_libraries(MBO gcov)
  endif()
endif()

install(TARGETS MBO DESTINATION lib)
