#include <stdlib.h>

#include "TensorOp.h"

/**
 * @brief Elementary embedding of an operator into a product space.
 * */
struct Embedding
{
	int i;
	ElemOp op;
};

/**
 * @brief A simple product of embeddings
 *
 * A SimpleTOp represents a single term in a sum making up a TOp.
 * */
struct SimpleTOp
{
	int numFactors;
	struct Embedding *embeddings;
};

void createEmbedding(struct Embedding *e);
void destroyEmbedding(struct Embedding *e);
static int findEmbedding(int i, int numEmbeddings,
			 struct Embedding *embeddings);
static void copyEmbedding(struct Embedding *dest, struct Embedding *src);
/*
 * @brief Search for an embedding.
 *
 * Gathers all embeddings embedded at space i into on slot and returns the index
 * of this slot.  This function reorders entries in embeddings.  Pointers into
 * this array are thus invalidated.
 *
 * @param i The slot to search for.
 * @param numEbeddings Number of embeddings in array.
 * @param embeddings Array of embeddings.
 * @return Index into array the embeddings array where resulting embeddings can
 * be found or a negative number if no suitable embedding has been found.
 * */
static int gatherIthEmbedding(int i, int numEmbeddings,
			      struct Embedding **embeddings);
static void destroySimpleTOp(struct SimpleTOp *term);
static void multiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);
static void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);

/**
   Data structure for tensor product operators.
   */
struct TensorOp
{
	/* Product space on which the operator is defined. */
	ProdSpace space;
	/* Number of terms (SimpleTOp's) making up the operator */
	int numTerms;
	/* Array of terms in sum */
	struct SimpleTOp *sum;
};

void tensorOpCreate(ProdSpace h, TensorOp *op)
{
	TensorOp a = malloc(sizeof(*a));
	a->space = prodSpaceCopy(h);
	a->numTerms = 0;
	a->sum = 0;
	*op = a;
}

void tensorOpDestroy(TensorOp *op)
{
	int i;
	for (i = 0; i < (*op)->numTerms; ++i) {
		destroySimpleTOp(&(*op)->sum[i]);
	}
	free((*op)->sum);
	prodSpaceDestroy(&(*op)->space);
	free(*op);
	*op = 0;
}

void destroySimpleTOp(struct SimpleTOp *term)
{
	int i;
	for (i = 0; i < term->numFactors; ++i) {
		destroyEmbedding(term->embeddings + i);
	}
	free(term->embeddings);
	free(term);
}

void tensorOpAddTo(ElemOp a, int i, TensorOp op)
{
	struct SimpleTOp *newTerm;
	struct Embedding *aEmbedded;

	op->sum = realloc(op->sum, (op->numTerms + 1) * sizeof(*op->sum));
	newTerm = op->sum;
	newTerm->numFactors = 1;
	newTerm->embeddings = malloc(sizeof(struct Embedding));
	aEmbedded = newTerm->embeddings;
	aEmbedded->i = i;
	aEmbedded->op = elemOpCopy(a);
	++op->numTerms;
}

void tensorOpAddScaledTo(double alpha, ElemOp a, int i, TensorOp op)
{
	struct SimpleTOp *newTerm;
	struct Embedding *aEmbedded;

	op->sum = realloc(op->sum, (op->numTerms + 1) * sizeof(*op->sum));
	newTerm = op->sum;
	newTerm->numFactors = 1;
	newTerm->embeddings = malloc(sizeof(struct Embedding));
	aEmbedded = newTerm->embeddings;
	aEmbedded->i = i;
	aEmbedded->op = elemOpCopy(a);
	elemOpScale(alpha, aEmbedded->op);
	++op->numTerms;
}

void tensorOpMul(TensorOp a, TensorOp *b)
{
	int i, j;
	int N = prodSpaceSize((*b)->space);

	for (i = 0; i < a->numTerms; ++i) {
		for (j = 0; j < (*b)->numTerms; ++j) {
			multiplySimpleTOps(N, a->sum + i, (*b)->sum + j);
		}
	}
}

void tensorOpPlus(TensorOp a, TensorOp *b)
{
	int i;
	(*b)->sum = realloc((*b)->sum, ((*b)->numTerms + a->numTerms) *
					   sizeof(*(*b)->sum));
	for (i = 0; i < a->numTerms; ++i) {
		copySimpleTOp((*b)->sum + (*b)->numTerms + i, a->sum + i);
	}
	(*b)->numTerms += a->numTerms;
}

void multiplySimpleTOps(int N, struct SimpleTOp *a, struct SimpleTOp *b)
{
	int i, j, k;
	for (i = 0; i < N; ++i) {
		j = gatherIthEmbedding(i, a->numFactors, &a->embeddings);
		k = gatherIthEmbedding(i, b->numFactors, &b->embeddings);
		if (j < 0) {
			/* Identity in a at i, needn't do anything */
		} else {
			if (k < 0) {
				/* Identity in b at i, but non-trivial operator
				   in a */
				b->embeddings =
				    realloc(b->embeddings,
					    (b->numFactors + 1) *
						sizeof(struct Embedding));
				copyEmbedding(b->embeddings + b->numFactors,
						a->embeddings + j);
				++b->numFactors;
			} else {
				/* Non-trivial operators in both slots; just
				   multiply them together. */
				elemOpMul(a->embeddings[j].op,
						&(b->embeddings[k].op));
			}
		}
	}
}

void copySimpleTOp(struct SimpleTOp* dest, struct SimpleTOp *src)
{
	int i;
	dest->numFactors = src->numFactors;
	dest->embeddings = realloc(
	    dest->embeddings, dest->numFactors * sizeof(*dest->embeddings));
	for (i = 0; i < src->numFactors; ++i) {
		copyEmbedding(dest->embeddings + i, src->embeddings + i);
	}
}

void copyEmbedding(struct Embedding *dest, struct Embedding *src)
{
	dest->i = src->i;
	elemOpDestroy(&dest->op);
	dest->op = elemOpCopy(src->op);
}

void destroyEmbedding(struct Embedding *e)
{
	elemOpDestroy(&e->op);
	free(e);
}

int gatherIthEmbedding(int i, int numEmbeddings, struct Embedding **embeddings)
{
	int first, current, next, numRemaining, j, jt;
	struct Embedding *tmp;
	first = findEmbedding(i, numEmbeddings, *embeddings);
	if (first < 0) return first;
	/* Accumulate all remaining Embeddings with the same slot into the first
	   one. Mark each further occurrance for deletion by setting its slot i
	   to a negative number. */
	current = first + 1;
	do {
		next = findEmbedding(i, numEmbeddings - current,
				     *embeddings + current);
		if (next > 0) {
			elemOpPlus((*embeddings)[next].op,
				   &(*embeddings)[first].op);
			(*embeddings)[next].i = -1;
		}
		current = next + 1;
	} while (next > 0);
	/* Count unmarked embeddings. */
	numRemaining = 0;
	for (j = 0; j < numEmbeddings; ++j) {
		if ((*embeddings)[j].i >= 0) {
			++numRemaining;
		}
	}
	/* Copy unmarked Embeddings into a new array. */
	tmp = malloc(numRemaining * sizeof(*tmp));
	jt = 0;
	for (j = 0; j < numEmbeddings; ++j) {
		if ((*embeddings)[j].i >= 0) {
			copyEmbedding(tmp + jt, (*embeddings) + j);
			++jt;
		}
	}
	/* Swap embeddings array with new array. */
	for (j = 0; j < numEmbeddings; ++j) {
		destroyEmbedding((*embeddings) + j);
	}
	free(*embeddings);
	*embeddings = tmp;
	return first;
}

int findEmbedding(int i, int numEmbeddings, struct Embedding *emb)
{
	int j;
	for (j = 0; j < numEmbeddings; ++j) {
		if (emb[j].i == i) return j;
	}
	return -1;
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
//
//static int testTensorOpCreate()
//{
//	int errs = 0;
//	TensorOp op;
//	tensorOpCreate(0, &op);
//	CHK_EQUAL(op->space, 0, errs);
//	CHK_EQUAL(op->sum, 0, errs);
//	tensorOpDestroy(&op);
//	return errs;
//}
//
//static int testTensorOpAddTo()
//{
//	int errs = 0;
//	TensorOp op;
//	ElemOp eop;
//	ProdSpace h;
//	double matrixElement = -3.4;
//
//	h = prodSpaceCreate(20);
//
//	elemOpCreate(&eop);
//	elemOpAddTo(14, 15, matrixElement, &eop);
//
//	tensorOpCreate(0, &op);
//	tensorOpAddTo(eop, 0, op);
//	CHK_EQUAL(op->sum->embedding->i, 0, errs);
//	CHK_EQUAL(elemOpCheck(op->sum->embedding->op), 0, errs);
//	CHK_EQUAL(op->sum->next, 0, errs);
//
//	tensorOpDestroy(&op);
//	elemOpDestroy(&eop);
//	prodSpaceDestroy(&h);
//	return errs;
//}
//
//static int testTensorOpAddScaledTo()
//{
//	int errs = 0;
//	TensorOp op;
//	ElemOp eop;
//	ProdSpace h;
//	double alpha = 2.1;
//
//	h = prodSpaceCreate(20);
//
//	elemOpCreate(&eop);
//	elemOpAddTo(14, 15, -3.4, &eop);
//
//	tensorOpCreate(0, &op);
//	tensorOpAddScaledTo(alpha, eop, 0, op);
//	CHK_EQUAL(op->sum->embedding->i, 0, errs);
//	CHK_EQUAL(elemOpCheck(op->sum->embedding->op), 0, errs);
//	CHK_EQUAL(op->sum->next, 0, errs);
//
//	tensorOpDestroy(&op);
//	elemOpDestroy(&eop);
//	prodSpaceDestroy(h);
//	return errs;
//}
//
//static int testFindEmbedding()
//{
//	int errs = 0;
//	int i;
//	int N = 10;
//	TensorOp op;
//	ElemOp eop;
//	ProdSpace h;
//	struct Embedding *emb;
//
//	h = prodSpaceCreate(20);
//	for (i = 0; i < N; ++i) {
//		prodSpaceMul(h, &h);
//	}
//
//	elemOpCreate(&eop);
//	elemOpAddTo(14, 15, -3.4, &eop);
//	elemOpAddTo(1, 5, 2.0, &eop);
//	elemOpAddTo(2, 3, 2.0, &eop);
//
//	tensorOpCreate(h, &op);
//	tensorOpAddTo(eop, 0, op);
//	tensorOpAddTo(eop, 3, op);
//
//	emb = findEmbedding(0, op->sum->next->embedding);
//	CHK_NOT_EQUAL(emb, 0, errs);
//	emb = findEmbedding(1, op->sum->next->embedding);
//	CHK_EQUAL(emb, 0, errs);
//	emb = findEmbedding(3, op->sum->embedding);
//	CHK_NOT_EQUAL(emb, 0, errs);
//	emb = findEmbedding(1, op->sum->embedding);
//	CHK_EQUAL(emb, 0, errs);
//
//	tensorOpDestroy(&op);
//	elemOpDestroy(&eop);
//	prodSpaceDestroy(&h);
//
//	return errs;
//}
//
//static int testTensorOpMul()
//{
//	int errs = 0;
//	int i;
//	int N = 10;
//	TensorOp op1, op2;
//	ElemOp eop1, eop2;
//	ProdSpace h;
//
//	h = prodSpaceCreate(20);
//	for (i = 0; i < N; ++i) {
//		prodSpaceMul(h, &h);
//	}
//
//	elemOpCreate(&eop1);
//	elemOpAddTo(14, 15, -3.4, &eop1);
//	elemOpAddTo(1, 5, 2.0, &eop1);
//	elemOpAddTo(2, 3, 2.0, &eop1);
//
//	elemOpCreate(&eop2);
//	elemOpAddTo(8, 3, -3.4, &eop2);
//	elemOpAddTo(3, 4, 2.0, &eop2);
//
//	tensorOpCreate(h, &op1);
//	tensorOpAddTo(eop1, 0, op1);
//	tensorOpCreate(h, &op2);
//	tensorOpAddTo(eop2, 0, op2);
//
//	tensorOpMul(op1, &op2);
//
//	CHK_EQUAL(op1->sum->embedding->i, 0, errs);
//	CHK_EQUAL(op1->sum->next, 0, errs);
//	CHK_EQUAL(elemOpCheck(op1->sum->embedding->op), 0, errs);
//
//	CHK_EQUAL(op2->sum->embedding->i, 0, errs);
//	CHK_EQUAL(op2->sum->embedding->next, 0, errs);
//	CHK_EQUAL(elemOpCheck(op2->sum->embedding->op), 0, errs);
//
//	tensorOpDestroy(&op1);
//	tensorOpDestroy(&op2);
//	elemOpDestroy(&eop1);
//	elemOpDestroy(&eop2);
//	prodSpaceDestroy(&h);
//	return errs;
//}
//
//static int testTensorOpPlus()
//{
//	int i, errs = 0, N = 5;
//	ProdSpace hTot, h1;
//	TensorOp op1, op2;
//	ElemOp sz, sp, sm;
//
//	hTot = prodSpaceCreate(0);
//	h1 = prodSpaceCreate(2);
//	for (i = 0; i < N; ++i) {
//		prodSpaceMul(h1, &hTot);
//	}
//
//	tensorOpCreate(hTot, &op1);
//	tensorOpCreate(hTot, &op2);
//	sz = sigmaZ();
//	sp = sigmaPlus();
//	sm = sigmaMinus();
//	tensorOpAddTo(sz, 0, op1);
//	tensorOpAddTo(sz, 1, op1);
//	tensorOpAddTo(sp, 0, op2);
//	tensorOpAddTo(sp, 2, op2);
//	tensorOpAddTo(sp, 3, op2);
//	tensorOpAddTo(sm, 3, op2);
//	tensorOpAddTo(sm, 2, op2);
//	CHK_EQUAL(0, tensorOpCheck(op1), errs);
//	CHK_EQUAL(0, tensorOpCheck(op2), errs);
//
//	tensorOpPlus(op1, &op2);
//	CHK_EQUAL(0, tensorOpCheck(op1), errs);
//
////	tensorOpPlus(op2, &op1);
////	CHK_EQUAL(0, tensorOpCheck(op1), errs);
//
//	elemOpDestroy(&sz);
//	elemOpDestroy(&sp);
//	elemOpDestroy(&sm);
//	tensorOpDestroy(&op1);
//	tensorOpDestroy(&op2);
//	prodSpaceDestroy(&hTot);
//	prodSpaceDestroy(&h1);
//	return errs;
//}
//
int tensorOpTest()
{
	int errs = 0;
//	errs += testTensorOpCreate();
//	errs += testTensorOpAddTo();
//	errs += testTensorOpAddScaledTo();
//	errs += testFindEmbedding();
//	errs += testTensorOpMul();
//	errs += testTensorOpPlus();
	return errs;
}
