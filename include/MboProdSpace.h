/**
 * @file MboProdSpace.h
 * @brief The MboProdSpace API.
 *
 * MboProdSpaces are used to describe tensor product spaces.
 * */
#ifndef MBO_PROD_SPACE_H
#define MBO_PROD_SPACE_H

#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboProdSpace;
/**
 * @brief Type for describing tensor product spaces
 * */
typedef struct MboProdSpace *MboProdSpace;

MboProdSpace mboProdSpaceCreate(MboLocInd);
void mboProdSpaceDestroy(MboProdSpace *);
void mboProdSpaceMul(MboProdSpace, MboProdSpace *);
MboProdSpace mboProdSpaceCopy(MboProdSpace);
MboGlobInd mboProdSpaceDim(MboProdSpace);
int mboProdSpaceSize(MboProdSpace);
void mboProdSpaceGetDims(MboProdSpace, int, MboLocInd *);
int mboProdSpaceEqual(MboProdSpace, MboProdSpace);
int mboProdSpaceCheck(MboProdSpace);
int mboProdSpaceTest();

#ifdef __cplusplus
}
#endif
#endif
