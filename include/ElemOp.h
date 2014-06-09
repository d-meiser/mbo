#ifndef ELEM_OP_H
#define ELEM_OP_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElemOp;
typedef struct ElemOp *ElemOp;

void CreateElemOp(ElemOp *);
void DestroyElemOp(ElemOp *);
void AddToElemOp(int, int, double, ElemOp *);
void ScaleElemOp(double, ElemOp);
void PlusElemOp(ElemOp, ElemOp *);
void MulElemOp(ElemOp, ElemOp *);

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
