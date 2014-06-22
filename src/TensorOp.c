#include <stdlib.h>

#include "TensorOp.h"

struct Embedding
{
	ElemOp op;
	int i;
	struct Embedding *next;
};

static struct Embedding *FindEmbedding(int i, struct Embedding *list);

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

static struct Embedding *GatherIthEmbedding(int, struct SimpleTOp *);
static void DestroySimpleTOp(struct SimpleTOp *term);
static void MultiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);

/**
   Data structure for tensor product operators.
   */
struct TensorOp
{
	/* Product space on which the operator is defined. */
	ProdSpace space;
	/* Pointer to first SimpleTOp in list making up the sum over SimpleTOp
	   operators. */
	struct SimpleTOp *sum;
};


void tensorOpCreate(ProdSpace h, TensorOp *op)
{
	TensorOp a = malloc(sizeof(*a));
	a->space = h;
	a->sum = 0;
	*op = a;
}

void tensorOpDestroy(TensorOp *op)
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

void tensorOpAddTo(ElemOp a, int i, TensorOp op)
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

void tensorOpAddScaledTo(double alpha, ElemOp a, int i, TensorOp op)
{
	tensorOpAddTo(a, i, op);
	ScaleElemOp(alpha, op->sum->embedding->op);
}

void tensorOpMul(TensorOp a, TensorOp *b)
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

int tensorOpCheck(TensorOp op)
{
	return 0;
}

/*
 * Tests
 * */
#include "TestUtils.h"
#include "ElemOp.h"

#define EPS 1.0e-12

int test_TensorOpCreate()
{
	int errs = 0;
	TensorOp op;
	tensorOpCreate(0, &op);
	CHK_EQUAL(op->space, 0, errs);
	CHK_EQUAL(op->sum, 0, errs);
	tensorOpDestroy(&op);
	return errs;
}

int test_TensorOpAddTo()
{
	int errs = 0;
	TensorOp op;
	ElemOp eop;
	ProdSpace h;
	double matrixElement = -3.4;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, matrixElement, &eop);

	tensorOpCreate(0, &op);
	tensorOpAddTo(eop, 0, op);
	CHK_EQUAL(op->sum->embedding->i, 0, errs);
	CHK_EQUAL(checkElemOp(op->sum->embedding->op), 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	tensorOpDestroy(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_TensorOpAddScaledTo()
{
	int errs = 0;
	TensorOp op;
	ElemOp eop;
	ProdSpace h;
	double alpha = 2.1;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, -3.4, &eop);

	tensorOpCreate(0, &op);
	tensorOpAddScaledTo(alpha, eop, 0, op);
	CHK_EQUAL(op->sum->embedding->i, 0, errs);
	CHK_EQUAL(checkElemOp(op->sum->embedding->op), 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	tensorOpDestroy(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_FindEmbedding()
{
	int errs = 0;
	int i;
	int N = 10;
	TensorOp op;
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

	tensorOpCreate(h, &op);
	tensorOpAddTo(eop, 0, op);
	tensorOpAddTo(eop, 3, op);

	emb = FindEmbedding(0, op->sum->next->embedding);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->next->embedding);
	CHK_EQUAL(emb, 0, errs);
	emb = FindEmbedding(3, op->sum->embedding);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->embedding);
	CHK_EQUAL(emb, 0, errs);

	tensorOpDestroy(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);

	return errs;
}

int test_TensorOpMul()
{
	int errs = 0;
	int i;
	int N = 10;
	TensorOp op1, op2;
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

	tensorOpCreate(h, &op1);
	tensorOpAddTo(eop1, 0, op1);
	tensorOpCreate(h, &op2);
	tensorOpAddTo(eop2, 0, op2);

	tensorOpMul(op1, &op2);

	CHK_EQUAL(op1->sum->embedding->i, 0, errs);
	CHK_EQUAL(op1->sum->next, 0, errs);
	CHK_EQUAL(checkElemOp(op1->sum->embedding->op), 0, errs);

	CHK_EQUAL(op2->sum->embedding->i, 0, errs);
	CHK_EQUAL(op2->sum->embedding->next, 0, errs);
	CHK_EQUAL(checkElemOp(op2->sum->embedding->op), 0, errs);

	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	DestroyElemOp(&eop1);
	DestroyElemOp(&eop2);
	DestroyProdSpace(&h);
	return errs;
}

int tensorOpTest()
{
	int errs = 0;
	errs += test_TensorOpCreate();
	errs += test_TensorOpAddTo();
	errs += test_TensorOpAddScaledTo();
	errs += test_FindEmbedding();
	errs += test_TensorOpMul();
	return errs;
}
