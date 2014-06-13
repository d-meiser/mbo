#include "TensProdOpImpl.h"

#include <stdlib.h>

static void Destroy(struct Embedding *sum);

void CreateTOp(ProdSpace h, TOp *op)
{
	TOp a = malloc(sizeof(*a));
	a->space = h;
	a->sum = 0;
	*op = a;
}

void DestroyTOp(TOp *op)
{
	struct SimpleTOp *nextTerm;

	if (*op == 0) return;
	DestroyProdSpace(&(*op)->space);
	while ((*op)->sum) {
		nextTerm = (*op)->sum;
		DestroySimpleTOp((*op)->sum);
		(*op)->sum = nextTerm;
	}
	free(*op);
	*op = 0;
}

void DestroySimpleTOp(struct SimpleTOp *term)
{
	struct Embedding* nextEmbedding;

	if (term == 0) return;
	while (term->embedding) {
		nextEmbedding = term->embedding->next;
		DestroyElemOp(&term->embedding->op);
		free(term->embedding);
		term->embedding = nextEmbedding;
	}
}

void AddToTOp(ElemOp a, int i, TOp op)
{
	struct Embedding *aEmbedded;
	struct SimpleTOp *simpleTOp;

	aEmbedded = malloc(sizeof(*aEmbedded));
	aEmbedded->op = a;
	aEmbedded->i = i;
	aEmbedded->next = 0;
	simpleTOp = malloc(sizeof(*simpleTOp));
	simpleTOp->embedding = aEmbedded;
	simpleTOp->next = op->sum;
	op->sum = simpleTOp;
}

void AddScaledToTOp(double alpha, ElemOp a, int i, TOp op)
{
	AddToTOp(a, i, op);
	ScaleElemOp(alpha, op->sum->embedding->op);
}

void MulTOp(TOp a, TOp *b)
{
	struct SimpleTOp* sa;
	struct SimpleTOp* sb;
	int N = SizeProdSpace((*b)->space);

	for (sa = a->sum; sa != 0; sa = sa->next) {
		for (sb = (*b)->sum; sb != 0; sb = sb->next) {
			MultiplySimpleTOps(N, sa, sb);
		}
	}
}

void MultiplySimpleTOps(int N, struct SimpleTOp *a, struct SimpleTOp *b)
{
	int i;
	struct Embedding *ea;
	struct Embedding *eb;

	for (i = 0; i < N; ++i) {
		ea = GatherIthEmbedding(i, a);
		eb = GatherIthEmbedding(i, b);
		if (ea != 0) {
			MulElemOp(ea->op, &eb->op);
		} else {
			/* Identity in ea, needn't do anything */
		}
	}
}

struct Embedding *GatherIthEmbedding(int i, struct SimpleTOp *sum)
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

struct Embedding *FindEmbedding(int i, struct Embedding *emb)
{
	while (emb) {
		if (emb->i == i) return emb;
		emb = emb->next;
	}
	return 0;
}

