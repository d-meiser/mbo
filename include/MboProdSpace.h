/**
 * @file MboProdSpace.h
 * @brief The MboProdSpace API.
 *
 * MboProdSpaces are used to describe tensor product spaces.
 * */
#ifndef MBO_PROD_SPACE_H
#define MBO_PROD_SPACE_H

#include <MboSys.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboProdSpace;
/**
 * @brief Type for describing tensor product spaces
 * */
typedef struct MboProdSpace *MboProdSpace;

MBO_API MboProdSpace mboProdSpaceCreate(MboLocInd);
MBO_API void mboProdSpaceDestroy(MboProdSpace *);
MBO_API void mboProdSpaceMul(MboProdSpace, MboProdSpace *);
MBO_API MboProdSpace mboProdSpaceCopy(MboProdSpace);
MBO_API MboGlobInd mboProdSpaceDim(MboProdSpace);
MBO_API int mboProdSpaceSize(MboProdSpace);
MBO_API void mboProdSpaceGetDims(MboProdSpace, int, MboLocInd *);
MBO_API int mboProdSpaceEqual(MboProdSpace, MboProdSpace);
MBO_API int mboProdSpaceCheck(MboProdSpace);
MBO_API int mboProdSpaceTest();

#ifdef __cplusplus
}
#endif
#endif
