#include "ElemOpImpl.h"
#include <stdlib.h>
#include <math.h>

double sqrt(double);

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

void PlusElemOp(ElemOp a, ElemOp *b)
{
	while (a) {
		AddToElemOp(a->op.m, a->op.n, a->op.val, b);
		a = a->next;
	}
}

ElemOp sigmaPlus()
{
	ElemOp sp;
	CreateElemOp(&sp);
	AddToElemOp(1, 0, 1.0, &sp);
	return sp;
}

ElemOp sigmaMinus()
{
	ElemOp sm;
	CreateElemOp(&sm);
	AddToElemOp(0, 1, 1.0, &sm);
	return sm;
}

ElemOp sigmaZ()
{
	ElemOp sz;
	CreateElemOp(&sz);
	AddToElemOp(1, 1, 1.0, &sz);
	AddToElemOp(0, 0, -1.0, &sz);
	return sz;
}

ElemOp eye(int d)
{
	ElemOp e;
        int i;
	CreateElemOp(&e);
	for (i = 0; i < d; ++i) {
		AddToElemOp(i, i, 1.0, &e);
	}
	return e;
}

ElemOp numOp(int d)
{
	ElemOp n;
        int i;
	CreateElemOp(&n);
	for (i = 1; i < d; ++i) {
		AddToElemOp(i, i, i, &n);
	}
	return n;
}

ElemOp annihilationOp(int d)
{
	ElemOp a;
        int i;
	CreateElemOp(&a);
	for (i = 1; i < d; ++i) {
		AddToElemOp(i - 1, i, sqrt(i), &a);
	}
	return a;
}

ElemOp creationOp(int d)
{
	ElemOp ad;
        int i;
	CreateElemOp(&ad);
	for (i = 1; i < d; ++i) {
		AddToElemOp(i, i - 1, sqrt(i), &ad);
	}
	return ad;
}


