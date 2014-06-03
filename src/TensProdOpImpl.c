#include "TensProdOpImpl.h"

#include <stdlib.h>

TOp CreateTOp(ProdSpace *h)
{
  TOp op = (TOp)malloc(sizeof(*op));
  return op;
}

void DestroyTOp(TOp op)
{
  free(op);
}

void AddToOp(ElemOp a, int i, TOp op)
{
  struct Embedding *next = op->sum;
  struct Embedding *aprime = malloc(sizeof(*aprime));
  aprime->op = a;
  aprime->i = i;
  aprime->next = next;
  op->sum = aprime;
}
