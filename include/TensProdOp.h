#ifndef TENS_PROD_OP_H
#define TENS_PROD_OP_H

#include "ProdSpace.h"
#include "ElemOp.h"

struct TOp;
typedef struct TOp *TOp;

void CreateTOp(ProdSpace, TOp *);
void DestroyTOp(TOp *);

void AddToTOp(ElemOp, int, TOp);

#endif
