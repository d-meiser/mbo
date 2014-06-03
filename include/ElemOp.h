#ifndef ELEM_OP_H
#define ELEM_OP_H

#include "BraKet.h"

struct ElemOp;
typedef struct ElemOp *ElemOp;

void CreateElemOp(struct BraKet, ElemOp *);
void DestroyElemOp(ElemOp);
void AddToElemOp(ElemOp, ElemOp);

#endif
