#include "TensProdOpImpl.h"

#include <stdlib.h>

static void Destroy(struct Embedding *head);

void CreateTOp(ProdSpace h, TOp *op)
{
	TOp a = malloc(sizeof(*a));
        a->space = h;
        a->sum = 0;
        *op = a;
}

void DestroyTOp(TOp *op)
{
	if (*op) {
		Destroy((*op)->sum);
		free(*op);
		*op = 0;
	}
}

void Destroy(struct Embedding *head)
{
	if (head) {
		Destroy(head->next);
		free(head);
	}
}

void AddToTOp(ElemOp a, int i, TOp op)
{
	struct Embedding *aprime = malloc(sizeof(*aprime));
	aprime->op = a;
	aprime->i = i;
	aprime->next = op->sum;
	op->sum = aprime;
}

void AddScaledToTOp(double alpha, ElemOp a, int i, TOp op)
{
	AddToTOp(a, i, op);
	ScaleElemOp(alpha, op->sum->op);
}

