#ifndef TENS_OP_H
#define TENS_OP_H

#include "ProdSpace.h"
#include "ElemOp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TensorOp;
typedef struct TensorOp *TensorOp;

void tensorOpCreate(ProdSpace, TensorOp *);
void tensorOpDestroy(TensorOp *);
void tensorOpAddTo(ElemOp, int, TensorOp);
void tensorOpAddScaledTo(double, ElemOp, int, TensorOp);

void tensorOpMul(TensorOp, TensorOp *);
void tensorOpPlus(TensorOp, TensorOp *);
void tensorOpScale(double alpha, TensorOp *);
void tensorOpKron(TensorOp, TensorOp, TensorOp *);
int tensorOpCheck(TensorOp);
int tensorOpTest();

#ifdef __cplusplus
}
#endif
#endif
