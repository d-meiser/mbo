set(QUO_HEADERS
    Mbo.h
    MboAmplitude.h
    MboElemOp.h
    MboErrors.h
    MboNonZeroEntry.h
    MboNumOp.h
    MboNumSubMatrix.h
    MboProdSpace.h
    MboTensorOp.h
    MboVec.h
   )

install(FILES ${QUO_HEADERS} ${PROJECT_BINARY_DIR}/MboExport.h DESTINATION include)

