#ifndef ELEM_OP_H
#define ELEM_OP_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElemOp;
typedef struct ElemOp *ElemOp;

void CreateElemOp(ElemOp *);
void DestroyElemOp(ElemOp *);
void AddToElemOp(int m, int n, double val, ElemOp *);

#ifdef __cplusplus
}
#endif

#endif
