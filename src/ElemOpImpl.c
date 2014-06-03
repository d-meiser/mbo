#include "ElemOpImpl.h"
#include <stdlib.h>

void CreateElemOp(struct BraKet bk, ElemOp *op)
{
	ElemOp newop = malloc(sizeof(*newop));
	newop->op = bk;
	newop->next = 0;
	*op = newop;
}

void DestroyElemOp(ElemOp op)
{
	if (op->next) {
		DestroyElemOp(op->next);
	}
	free(op);
}

void AddToElemOp(ElemOp a, ElemOp b)
{
	ElemOp last = ElemOpLast(b);
	ElemOp aprime = ElemOpCopy(a);
	last->next = aprime;
}

ElemOp ElemOpLast(ElemOp op)
{
	ElemOp last = op;
	while (last->next) {
		last = last->next;
	}
	return last;
}

ElemOp ElemOpCopy(ElemOp a)
{
	ElemOp copy = 0;
	if (a) {
		copy = malloc(sizeof(*copy));
		copy->op = a->op;
		copy->next = ElemOpCopy(a->next);
	}
	return copy;
}
