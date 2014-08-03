#ifndef MBO_PROD_SPACE_H
#define MBO_PROD_SPACE_H

#ifdef __cplusplus
extern "C" {
#endif

struct MboProdSpace;
typedef struct MboProdSpace *MboProdSpace;

MboProdSpace mboProdSpaceCreate(int);
void mboProdSpaceDestroy(MboProdSpace *);
void mboProdSpaceMul(MboProdSpace, MboProdSpace *);
MboProdSpace mboProdSpaceCopy(MboProdSpace);
long long mboProdSpaceDim(MboProdSpace);
int mboProdSpaceSize(MboProdSpace);
void mboProdSpaceGetDims(MboProdSpace, int, int *);
int mboProdSpaceEqual(MboProdSpace, MboProdSpace);
int mboProdSpaceCheck(MboProdSpace);
int mboProdSpaceTest();

#ifdef __cplusplus
}
#endif
#endif
