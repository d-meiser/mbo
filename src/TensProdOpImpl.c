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

void MulTOp(TOp a, TOp *b)
{
	int i;
	int N = SizeProdSpace((*b)->space);
	for (i = 0; i < N; ++i) {
		MultiplyIthEmbeddings(i, a->sum, (*b)->sum);
	}
}

void MultiplyIthEmbeddings(int i, struct Embedding *asum,
			   struct Embedding *bsum)
{
	struct Embedding *ea = GatherIthEmbedding(i, asum);
	struct Embedding *eb = GatherIthEmbedding(i, bsum);
	if (ea != 0 && eb != 0) {
	}
	if (ea == 0 && eb != 0) {
	}
	if (ea != 0 && eb == 0) {
	}
	if (ea == 0 && eb == 0) {
	}
}

struct Embedding *FindEmbedding(int i, struct Embedding *sum)
{
	while (sum) {
		if (sum->i == i) return sum;
		sum = sum->next;
	}
	return 0;
}

struct Embedding *GatherIthEmbedding(int i, struct Embedding *sum)
{
	struct Embedding *first, *prev, *next;
	first = FindEmbedding(i, sum);
	if (first == 0) return first;
	prev = first;
	next = prev->next;
	while (next) {
		if (next->i == i) {
			PlusElemOp(next->op, &first->op);
			prev->next = next->next;
			free(next);
		}
		prev = prev->next;
		next = prev->next;
	}
	return first;
}

