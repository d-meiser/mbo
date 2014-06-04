#ifndef ELEM_OP_H
#define ELEM_OP_H

struct ElemOp;
typedef struct ElemOp *ElemOp;

void CreateElemOp(ElemOp *);
void DestroyElemOp(ElemOp *);
void AddToElemOp(int m, int n, double val, ElemOp *);

#endif
