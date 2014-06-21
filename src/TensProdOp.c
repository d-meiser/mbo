#include <stdlib.h>

#include "TensProdOp.h"
struct Embedding
{
	ElemOp op;
	int i;
	struct Embedding *next;
};

struct Embedding *FindEmbedding(int i, struct Embedding *list);
void MultiplyEmbeddings(int, struct Embedding *a, struct Embedding *b);

/**
 * @brief A simple product of embeddings
 *
 * A SimpleTOp represents a single term in a sum making up a TOp.
 * */
struct SimpleTOp
{
	struct Embedding *embedding;
	struct SimpleTOp *next;
};

struct Embedding *GatherIthEmbedding(int, struct SimpleTOp *);
void DestroySimpleTOp(struct SimpleTOp *term);
void MultiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);

/**
   Data structure for tensor product operators.
   */
struct TOp
{
	/* Product space on which the operator is defined. */
	ProdSpace space;
	/* Pointer to first SimpleTOp in list making up the sum over SimpleTOp
	   operators. */
	struct SimpleTOp *sum;
};


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
	while ((*op)->sum) {
		nextTerm = (*op)->sum->next;
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
	free(term);
}

void AddToTOp(ElemOp a, int i, TOp op)
{
	struct Embedding *aEmbedded;
	struct SimpleTOp *simpleTOp;

	aEmbedded = malloc(sizeof(*aEmbedded));
	CreateElemOp(&aEmbedded->op);
	PlusElemOp(a, &aEmbedded->op);
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
	first = FindEmbedding(i, sum->embedding);
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

int tensProdOpCheck(TOp op)
{
	return 0;
}

/*
 * Tests
 * */
#include "TestUtils.h"
#include "ElemOp.h"

#define EPS 1.0e-12

int test_CreateTOp()
{
	int errs = 0;
	TOp op;
	CreateTOp(0, &op);
	CHK_EQUAL(op->space, 0, errs);
	CHK_EQUAL(op->sum, 0, errs);
	DestroyTOp(&op);
	return errs;
}

int test_AddToTOp()
{
	int errs = 0;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	double matrixElement = -3.4;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, matrixElement, &eop);

	CreateTOp(0, &op);
	AddToTOp(eop, 0, op);
	CHK_EQUAL(op->sum->embedding->i, 0, errs);
	CHK_EQUAL(checkElemOp(op->sum->embedding->op), 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_AddScaledToTOp()
{
	int errs = 0;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	double alpha = 2.1;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, -3.4, &eop);

	CreateTOp(0, &op);
	AddScaledToTOp(alpha, eop, 0, op);
	CHK_EQUAL(op->sum->embedding->i, 0, errs);
	CHK_EQUAL(checkElemOp(op->sum->embedding->op), 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_FindEmbedding()
{
	int errs = 0;
	int i;
	int N = 10;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	struct Embedding *emb;

	h = CreateProdSpace(20);
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &h);
	}

	CreateElemOp(&eop);
	AddToElemOp(14, 15, -3.4, &eop);
	AddToElemOp(1, 5, 2.0, &eop);
	AddToElemOp(2, 3, 2.0, &eop);

	CreateTOp(h, &op);
	AddToTOp(eop, 0, op);
	AddToTOp(eop, 3, op);

	emb = FindEmbedding(0, op->sum->next->embedding);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->next->embedding);
	CHK_EQUAL(emb, 0, errs);
	emb = FindEmbedding(3, op->sum->embedding);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->embedding);
	CHK_EQUAL(emb, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);

	return errs;
}

int test_MulTOp()
{
	int errs = 0;
	int i;
	int N = 10;
	TOp op1, op2;
	ElemOp eop1, eop2;
	ProdSpace h;

	h = CreateProdSpace(20);
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &h);
	}

	CreateElemOp(&eop1);
	AddToElemOp(14, 15, -3.4, &eop1);
	AddToElemOp(1, 5, 2.0, &eop1);
	AddToElemOp(2, 3, 2.0, &eop1);

	CreateElemOp(&eop2);
	AddToElemOp(8, 3, -3.4, &eop2);
	AddToElemOp(3, 4, 2.0, &eop2);

	CreateTOp(h, &op1);
	AddToTOp(eop1, 0, op1);
	CreateTOp(h, &op2);
	AddToTOp(eop2, 0, op2);

	MulTOp(op1, &op2);

	CHK_EQUAL(op1->sum->embedding->i, 0, errs);
	CHK_EQUAL(op1->sum->next, 0, errs);
	CHK_EQUAL(checkElemOp(op1->sum->embedding->op), 0, errs);

	CHK_EQUAL(op2->sum->embedding->i, 0, errs);
	CHK_EQUAL(op2->sum->embedding->next, 0, errs);
	CHK_EQUAL(checkElemOp(op2->sum->embedding->op), 0, errs);

	DestroyTOp(&op1);
	DestroyTOp(&op2);
	DestroyElemOp(&eop1);
	DestroyElemOp(&eop2);
	DestroyProdSpace(&h);
	return errs;
}

int tensProdOpTest()
{
	int errs = 0;
	errs += test_CreateTOp();
	errs += test_AddToTOp();
	errs += test_AddScaledToTOp();
	errs += test_FindEmbedding();
	errs += test_MulTOp();
	return errs;
}
