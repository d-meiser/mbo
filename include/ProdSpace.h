#ifndef PROD_SPACE_H
#define PROD_SPACE_H

#ifdef __cplusplus
extern "C" {
#endif

struct ProdSpace;
typedef struct ProdSpace *ProdSpace;

ProdSpace prodSpaceCreate(int);
void prodSpaceDestroy(ProdSpace *);
void prodSpaceMul(ProdSpace, ProdSpace *);
ProdSpace prodSpaceCopy(ProdSpace);
long long prodSpaceDim(ProdSpace);
int prodSpaceSize(ProdSpace);
int prodSpaceEqual(ProdSpace, ProdSpace);
int prodSpaceCheck(ProdSpace);
int prodSpaceTest();

#ifdef __cplusplus
}
#endif
#endif
