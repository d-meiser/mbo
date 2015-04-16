include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    )

set(MBO_SRCS
  MboElemOp.c
  MboNumOp.c
  MboNumSubMatrix.c
  MboProdSpace.c
  MboTensorOp.c
  MboTensorOpPrivate.c
  MboVec.c
  Embedding.c
  Integrator.c
  SimpleTOp.c
  Tile.c
  Utilities.c
  )
set(MBO_HDRS
  Embedding.h
  Integrator.h
  MboNumOpPrivate.h
  MboTensorOpPrivate.h
  SimpleTOp.h
  Tile.h
  Utilities.h
  ../include/Mbo.h
  ../include/MboAmplitude.h
  ../include/MboElemOp.h
  ../include/MboErrors.h
  ../include/MboIndices.h
  ../include/MboNonZeroEntry.h
  ../include/MboNumOp.h
  ../include/MboNumSubMatrix.h
  ../include/MboProdSpace.h
  ../include/MboTensorOp.h
  ../include/MboVec.h
  )

add_compiler_export_flags()
add_library(MBO
    ${MBO_SRCS}
    ${MBO_HDRS}
    )
if(UNIX)
  if(NOT EXPLICIT_MATH_LIB_NOT_NEEDED)
    target_link_libraries(MBO m)
  endif()
  if(MBO_ENABLE_COVERAGE)
    target_link_libraries(MBO gcov)
  endif()
endif()

install(TARGETS MBO DESTINATION lib)
