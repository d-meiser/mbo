#include "TensProdOpImpl.h"

#include <stdlib.h>

static void Destroy(struct EmbeddingList *sum);

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

void Destroy(struct EmbeddingList *sum)
{
	struct Embedding* e;
	if (sum) {
		Destroy(sum->next);
		e = sum->first->next;
		while (e) {
			free(sum->first);
			e = e->next;
		}
		free(sum);
	}
}

void AddToTOp(ElemOp a, int i, TOp op)
{
	struct EmbeddingList *sum = malloc(sizeof(*sum));
	sum->first = 0;
	sum->next = op->sum;
	struct Embedding *aprime = malloc(sizeof(*aprime));
	aprime->op = a;
	aprime->i = i;
	aprime->next = 0;
	op->sum = sum;
}

void AddScaledToTOp(double alpha, ElemOp a, int i, TOp op)
{
	AddToTOp(a, i, op);
	ScaleElemOp(alpha, op->sum->first->op);
}

void MulTOp(TOp a, TOp *b)
{
	struct Embedding* ea;
	struct Embedding* eb;
	int N = SizeProdSpace((*b)->space);

	for (ea = a->sum; ea != 0; ea = ea->next) {
		for (eb = (*b)->sum; eb != 0; eb = eb->next) {
			MultiplyEmbeddings(N, ea, eb);
		}
	}
}

void MultiplyEmbeddings(int N, struct Embedding *asum, struct Embedding *bsum)
{
	int i;
	struct Embedding *ea;
	struct Embedding *eb;

	for (i = 0; i < N; ++i) {
		ea = GatherIthEmbedding(i, asum);
		eb = GatherIthEmbedding(i, bsum);
		if (ea != 0) {
			MulElemOp(ea->op, &eb->op);
		} else {
			/* Identity in ea, needn't do anything */
		}
	}
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

struct Embedding *FindEmbedding(int i, struct Embedding *emb)
{
	while (emb) {
		if (emb->i == i) return emb;
		emb = emb->next;
	}
	return 0;
}

