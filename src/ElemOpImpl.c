#include "ElemOpImpl.h"
#include <stdlib.h>

void CreateElemOp(ElemOp *op)
{
	*op = 0;
}

void DestroyElemOp(ElemOp *op)
{
	if (*op) {
		DestroyElemOp(&(*op)->next);
		free(*op);
		*op = 0;
	}
}

void AddToElemOp(int m, int n, double val, ElemOp *a)
{
	struct ElemOp *b = malloc(sizeof(*b));
	b->op.m = m;
	b->op.n = n;
	b->op.val = val;
	b->next = *a;
	*a = b;
}

void ScaleElemOp(double alpha, ElemOp op)
{
	ElemOp a;
	for (a = op; a != 0; a = a->next) {
		a->op.val *= alpha;
	}
}

