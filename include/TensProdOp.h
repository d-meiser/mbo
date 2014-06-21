#ifndef TENS_PROD_OP_H
#define TENS_PROD_OP_H

#include "ProdSpace.h"
#include "ElemOp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TOp;
typedef struct TOp *TOp;

void CreateTOp(ProdSpace, TOp *);
void DestroyTOp(TOp *);
void AddToTOp(ElemOp, int, TOp);
void AddScaledToTOp(double, ElemOp, int, TOp);

void MulTOp(TOp, TOp *);
void PlusOp(TOp, TOp *);
void ScaleTOp(double alpha, TOp *);
void KronTOp(TOp, TOp, TOp *);

#ifdef __cplusplus
}
#endif
#endif
