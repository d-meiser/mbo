/**
 * @file MboProdSpace.h
 * @brief The MboProdSpace API.
 *
 * MboProdSpaces are used to describe tensor product spaces.
 * */
#ifndef MBO_PROD_SPACE_H
#define MBO_PROD_SPACE_H

#include <MboExport.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboProdSpace_t;
/**
 * @brief Type for describing tensor product spaces
 * */
typedef struct MboProdSpace_t *MboProdSpace;

MBO_EXPORT MboProdSpace mboProdSpaceCreate(MboLocInd);
MBO_EXPORT void mboProdSpaceDestroy(MboProdSpace *);
MBO_EXPORT void mboProdSpaceMul(MboProdSpace, MboProdSpace *);
MBO_EXPORT MboProdSpace mboProdSpaceCopy(MboProdSpace);
MBO_EXPORT MboGlobInd mboProdSpaceDim(MboProdSpace);
MBO_EXPORT int mboProdSpaceSize(MboProdSpace);
MBO_EXPORT void mboProdSpaceGetDims(MboProdSpace, int, MboLocInd *);
MBO_EXPORT int mboProdSpaceEqual(MboProdSpace, MboProdSpace);
MBO_EXPORT int mboProdSpaceCheck(MboProdSpace);

#ifdef __cplusplus
}
#endif
#endif
