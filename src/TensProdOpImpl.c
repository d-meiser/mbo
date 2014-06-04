#include "TensProdOpImpl.h"

#include <stdlib.h>

void CreateTOp(ProdSpace h, TOp *op)
{
	TOp a = malloc(sizeof(*a));
        a->space = h;
        a->sum = 0;
        *op = a;
}

void DestroyTOp(TOp *op)
{
	free(*op);
        *op = 0;
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
