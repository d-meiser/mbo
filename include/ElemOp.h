#ifndef ELEM_OP_H
#define ELEM_OP_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElemOp;
typedef struct ElemOp *ElemOp;

void elemOpCreate(ElemOp *);
void elemOpDestroy(ElemOp *);
void elemOpAddTo(int, int, double, ElemOp *);
void elemOpScale(double, ElemOp);
void elemOpPlus(ElemOp, ElemOp *);
void elemOpMul(ElemOp, ElemOp *);
int elemOpCheck(ElemOp);
int elemOpTest();

ElemOp sigmaPlus();
ElemOp sigmaMinus();
ElemOp sigmaZ();
ElemOp eye(int);
ElemOp numOp(int);
ElemOp annihilationOp(int);
ElemOp creationOp(int);

#ifdef __cplusplus
}
#endif

#endif
