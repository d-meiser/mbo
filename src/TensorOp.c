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

static void destroyEmbedding(struct Embedding *e);
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
static int gatherIthEmbedding(int i, int *numEmbeddings,
			      struct Embedding **embeddings);
static void destroySimpleTOp(struct SimpleTOp *term);
static void multiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);
static void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);
static void scaleSimpleTOp(double alpha, ProdSpace h, struct SimpleTOp *op);
static int checkSimpleTOp(struct SimpleTOp *sa);

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

void tensorOpNull(ProdSpace h, TensorOp *op)
{
	TensorOp a = malloc(sizeof(*a));
	a->space = prodSpaceCopy(h);
	a->numTerms = 0;
	a->sum = 0;
	*op = a;
}

void tensorOpIdentity(ProdSpace h, TensorOp *op)
{
	TensorOp a = malloc(sizeof(*a));
	a->space = prodSpaceCopy(h);
	a->numTerms = 1;
	a->sum = malloc(sizeof(*a->sum));
	a->sum[0].numFactors = 0;
	a->sum[0].embeddings = 0;
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
}

void tensorOpAddTo(ElemOp a, int i, TensorOp op)
{
	struct SimpleTOp *newTerm;
	struct Embedding *aEmbedded;

	op->sum = realloc(op->sum, (op->numTerms + 1) * sizeof(*op->sum));
	newTerm = op->sum + op->numTerms;
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
	newTerm = op->sum + op->numTerms;
	newTerm->numFactors = 1;
	newTerm->embeddings = malloc(sizeof(struct Embedding));
	aEmbedded = newTerm->embeddings;
	aEmbedded->i = i;
	aEmbedded->op = elemOpCopy(a);
	elemOpScale(alpha, aEmbedded->op);
	++op->numTerms;
}

void tensorOpMul(TensorOp a, TensorOp b, TensorOp *c)
{
	int i, j;
	int N;
	int numNewTerms = a->numTerms * b->numTerms;
	struct SimpleTOp *cFillPtr;

	if (numNewTerms) {
		(*c)->sum = realloc((*c)->sum, ((*c)->numTerms + numNewTerms) *
						   (sizeof(*(*c)->sum)));
	}

	N = prodSpaceSize(b->space);
	cFillPtr = (*c)->sum + (*c)->numTerms;
	for (i = 0; i < a->numTerms; ++i) {
		for (j = 0; j < b->numTerms; ++j) {
			cFillPtr->numFactors = 0;
			cFillPtr->embeddings = 0;
			copySimpleTOp(cFillPtr, b->sum + j);
			multiplySimpleTOps(N, a->sum + i, cFillPtr);
			++cFillPtr;
		}
	}
	(*c)->numTerms += numNewTerms;
}

void tensorOpPlus(TensorOp a, TensorOp *b)
{
	int i;
	(*b)->sum = realloc((*b)->sum, ((*b)->numTerms + a->numTerms) *
					   sizeof(*(*b)->sum));
	for (i = 0; i < a->numTerms; ++i) {
		(*b)->sum[(*b)->numTerms + i].numFactors = 0;
		(*b)->sum[(*b)->numTerms + i].embeddings = 0;
		copySimpleTOp((*b)->sum + (*b)->numTerms + i, a->sum + i);
	}
	(*b)->numTerms += a->numTerms;
}

void tensorOpScale(double alpha, TensorOp *a)
{
	int i;
	for (i = 0; i < (*a)->numTerms; ++i) {
		scaleSimpleTOp(alpha, (*a)->space, (*a)->sum + i);
	}
}

void multiplySimpleTOps(int N, struct SimpleTOp *a, struct SimpleTOp *b)
{
	int i, j, k;
	for (i = 0; i < N; ++i) {
		j = gatherIthEmbedding(i, &a->numFactors, &a->embeddings);
		k = gatherIthEmbedding(i, &b->numFactors, &b->embeddings);
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
				elemOpCreate(&b->embeddings[b->numFactors].op);
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
		elemOpCreate(&dest->embeddings[i].op);
		copyEmbedding(dest->embeddings + i, src->embeddings + i);
	}
}
void scaleSimpleTOp(double alpha, ProdSpace h, struct SimpleTOp *op)
{
	int d;
	if (prodSpaceSize(h) == 0) return;
	if (op->numFactors == 0) {
		op->embeddings = malloc(sizeof(*op->embeddings));
		op->embeddings[0].i = 0;
		prodSpaceGetDims(h, 1, &d);
		op->embeddings[0].op = eye(d);
		op->numFactors += 1;
	}
	elemOpScale(alpha, op->embeddings[0].op);
}

void copyEmbedding(struct Embedding *dest, struct Embedding *src)
{
	dest->i = src->i;
	if (dest->op) {
		elemOpDestroy(&dest->op);
	}
	dest->op = elemOpCopy(src->op);
}

void destroyEmbedding(struct Embedding *e)
{
	elemOpDestroy(&e->op);
}

int gatherIthEmbedding(int i, int *numEmbeddings, struct Embedding **embeddings)
{
	int first, current, next, numRemaining, j, jt;
	struct Embedding *tmp;
	first = findEmbedding(i, *numEmbeddings, *embeddings);
	if (first < 0) return first;
	/* Accumulate all remaining Embeddings with the same slot into the first
	   one. Mark each further occurrance for deletion by setting its slot i
	   to a negative number. */
	current = first + 1;
	do {
		next = findEmbedding(i, *numEmbeddings - current,
				     *embeddings + current);
		if (next > 0) {
			next += current;
			elemOpPlus((*embeddings)[next].op,
				   &(*embeddings)[first].op);
			(*embeddings)[next].i = -1;
		}
		current = next + 1;
	} while (next > 0);
	/* Count unmarked embeddings. */
	numRemaining = 0;
	for (j = 0; j < *numEmbeddings; ++j) {
		if ((*embeddings)[j].i >= 0) {
			++numRemaining;
		}
	}
	/* Copy unmarked Embeddings into a new array. */
	tmp = malloc(numRemaining * sizeof(*tmp));
	jt = 0;
	for (j = 0; j < *numEmbeddings; ++j) {
		if ((*embeddings)[j].i >= 0) {
			tmp[jt].op = 0;
			copyEmbedding(tmp + jt, (*embeddings) + j);
			++jt;
		}
	}
	/* Swap embeddings array with new array. */
	for (j = 0; j < *numEmbeddings; ++j) {
		destroyEmbedding((*embeddings) + j);
	}
	free(*embeddings);
	*numEmbeddings = numRemaining;
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
	int i, errs = 0;
	if (op->numTerms < 0) ++errs;
	errs += prodSpaceCheck(op->space);
	for (i = 0; i < op->numTerms; ++i) {
		errs += checkSimpleTOp(op->sum + i);
	}
	return errs;
}

int checkSimpleTOp(struct SimpleTOp *sa)
{
	int i, errs = 0;
	if (sa->numFactors < 0) ++errs;
	for (i = 0; i < sa->numFactors; ++i) {
		if (sa->embeddings[i].i < 0) ++errs;
		errs += elemOpCheck(sa->embeddings[i].op);
	}
	return errs;
}

/*
 * Tests
 * */
#include "TestUtils.h"

#define EPS 1.0e-12

static int testTensorOpNull()
{
	int errs = 0;
	TensorOp op;
	ProdSpace h = prodSpaceCreate(1);
	tensorOpNull(h, &op);
	CHK_TRUE(prodSpaceEqual(op->space, h), errs);
	CHK_EQUAL(op->numTerms, 0, errs);
	tensorOpDestroy(&op);
	prodSpaceDestroy(&h);
	return errs;
}

static int testTensorOpIdentity()
{
	int errs = 0;
	TensorOp op;
	ProdSpace h = prodSpaceCreate(1);
	tensorOpIdentity(h, &op);
	CHK_TRUE(prodSpaceEqual(op->space, h), errs);
	CHK_EQUAL(op->numTerms, 1, errs);
	tensorOpDestroy(&op);
	prodSpaceDestroy(&h);
	return errs;
}

static int testTensorOpAddTo()
{
	int errs = 0;
	TensorOp op;
	ElemOp eop;
	ProdSpace h;
	double matrixElement = -3.4;

	h = prodSpaceCreate(20);

	elemOpCreate(&eop);
	elemOpAddTo(14, 15, matrixElement, &eop);

	tensorOpNull(h, &op);
	tensorOpAddTo(eop, 0, op);
	CHK_EQUAL(op->numTerms, 1, errs);
	CHK_TRUE(prodSpaceEqual(op->space, h), errs);
	CHK_TRUE(op->sum != 0, errs);
	CHK_EQUAL(op->sum[0].numFactors, 1, errs);
	CHK_TRUE(op->sum[0].embeddings != 0, errs);
	CHK_EQUAL(op->sum[0].embeddings[0].i, 0, errs);

	tensorOpDestroy(&op);
	elemOpDestroy(&eop);
	prodSpaceDestroy(&h);
	return errs;
}

static int testTensorOpAddScaledTo()
{
	int errs = 0;
	TensorOp op;
	ElemOp eop;
	ProdSpace h;
	double alpha = 2.1;

	h = prodSpaceCreate(20);

	elemOpCreate(&eop);
	elemOpAddTo(14, 15, -3.4, &eop);

	tensorOpNull(h, &op);
	tensorOpAddScaledTo(alpha, eop, 0, op);
	CHK_EQUAL(op->numTerms, 1, errs);
	CHK_TRUE(prodSpaceEqual(op->space, h), errs);
	CHK_TRUE(op->sum != 0, errs);
	CHK_EQUAL(op->sum[0].numFactors, 1, errs);
	CHK_TRUE(op->sum[0].embeddings != 0, errs);
	CHK_EQUAL(op->sum[0].embeddings[0].i, 0, errs);

	tensorOpAddScaledTo(alpha, eop, 1, op);
	CHK_EQUAL(op->numTerms, 2, errs);

	tensorOpDestroy(&op);
	elemOpDestroy(&eop);
	prodSpaceDestroy(&h);
	return errs;
}

static int testFindEmbedding()
{
	int errs = 0;
	int n = 5;
	struct Embedding *embeddings = malloc(n * sizeof(*embeddings));

	embeddings[0].i = 3;
	embeddings[1].i = 4;
	embeddings[2].i = -1;
	embeddings[3].i = 3;
	embeddings[4].i = 0;

	CHK_EQUAL(findEmbedding(0, n, embeddings), 4, errs);
	CHK_EQUAL(findEmbedding(3, n, embeddings), 0, errs);
	CHK_EQUAL(findEmbedding(3, n - 1, embeddings + 1), 2, errs);
	CHK_TRUE(findEmbedding(2, n, embeddings) < 0, errs);
	CHK_EQUAL(findEmbedding(-1, n, embeddings), 2, errs);

	free(embeddings);

	return errs;
}

static int testGatherIthEmbedding()
{
	int errs = 0;
	int i, n = 5;
	struct Embedding *embeddings = malloc(n * sizeof(*embeddings));

	embeddings[0].i = 3;
	elemOpCreate(&embeddings[0].op);
	embeddings[1].i = 4;
	elemOpCreate(&embeddings[1].op);
	embeddings[2].i = 0;
	elemOpCreate(&embeddings[2].op);
	embeddings[3].i = 3;
	elemOpCreate(&embeddings[3].op);
	embeddings[4].i = 0;
	elemOpCreate(&embeddings[4].op);

	i = gatherIthEmbedding(0, &n, &embeddings);
	CHK_EQUAL(i, 2, errs);
	CHK_EQUAL(n, 4, errs);

	i = gatherIthEmbedding(4, &n, &embeddings);
	CHK_EQUAL(i, 1, errs);
	CHK_EQUAL(n, 4, errs);

	i = gatherIthEmbedding(3, &n, &embeddings);
	CHK_EQUAL(i, 0, errs);
	CHK_EQUAL(n, 3, errs);

	for (i = 0; i < n; ++i) {
		elemOpDestroy(&embeddings[i].op);
	}

	free(embeddings);

	return errs;
}

static int testTensorOpMul()
{
	int errs = 0;
	int i;
	int N = 10;
	TensorOp op1, op2, op3;
	ElemOp eop1, eop2;
	ProdSpace h1, h2;

	h1 = prodSpaceCreate(5);
	h2 = prodSpaceCreate(0);
	for (i = 0; i < N; ++i) {
		prodSpaceMul(h1, &h2);
	}

	elemOpCreate(&eop1);
	elemOpAddTo(14, 15, -3.4, &eop1);
	elemOpAddTo(1, 5, 2.0, &eop1);
	elemOpAddTo(2, 3, 2.0, &eop1);

	elemOpCreate(&eop2);
	elemOpAddTo(8, 3, -3.4, &eop2);
	elemOpAddTo(3, 4, 2.0, &eop2);

	tensorOpNull(h2, &op1);
	tensorOpNull(h2, &op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	tensorOpNull(h2, &op1);
	tensorOpAddTo(eop1, 0, op1);
	tensorOpNull(h2, &op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	tensorOpNull(h2, &op1);
	tensorOpNull(h2, &op2);
	tensorOpAddTo(eop1, 0, op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	tensorOpNull(h2, &op1);
	tensorOpAddTo(eop1, 0, op1);
	tensorOpNull(h2, &op2);
	tensorOpAddTo(eop2, 0, op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 1, errs);
	CHK_EQUAL(op3->sum[0].numFactors, 1, errs);
	CHK_EQUAL(op3->sum[0].embeddings[0].i, 0, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	tensorOpNull(h2, &op1);
	tensorOpAddTo(eop1, 0, op1);
	tensorOpNull(h2, &op2);
	tensorOpAddTo(eop2, 1, op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 1, errs);
	CHK_EQUAL(op3->sum[0].numFactors, 2, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	tensorOpNull(h2, &op1);
	tensorOpAddTo(eop1, 0, op1);
	tensorOpAddTo(eop1, 1, op1);
	tensorOpNull(h2, &op2);
	tensorOpAddTo(eop2, 1, op2);
	tensorOpAddTo(eop2, 2, op2);
	tensorOpNull(h2, &op3);
	tensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 4, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	tensorOpDestroy(&op3);

	elemOpDestroy(&eop1);
	elemOpDestroy(&eop2);
	prodSpaceDestroy(&h1);
	prodSpaceDestroy(&h2);
	return errs;
}

static int testTensorOpPlus()
{
	int i, errs = 0, N = 5;
	ProdSpace hTot, h1;
	TensorOp op1, op2;
	ElemOp sz, sp, sm;

	hTot = prodSpaceCreate(0);
	h1 = prodSpaceCreate(2);
	for (i = 0; i < N; ++i) {
		prodSpaceMul(h1, &hTot);
	}

	sz = sigmaZ();
	sp = sigmaPlus();
	sm = sigmaMinus();

	tensorOpNull(hTot, &op1);
	tensorOpNull(hTot, &op2);
	tensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 0, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);

	tensorOpNull(hTot, &op1);
	tensorOpNull(hTot, &op2);
	tensorOpAddTo(sz, 1, op2);
	tensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 1, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);

	tensorOpNull(hTot, &op1);
	tensorOpAddTo(sz, 1, op1);
	tensorOpNull(hTot, &op2);
	tensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 1, errs);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);

	tensorOpNull(hTot, &op1);
	tensorOpNull(hTot, &op2);
	tensorOpAddTo(sz, 0, op1);
	tensorOpAddTo(sz, 1, op1);
	tensorOpAddTo(sp, 0, op2);
	tensorOpAddTo(sp, 2, op2);
	tensorOpAddTo(sp, 3, op2);
	tensorOpAddTo(sm, 3, op2);
	tensorOpAddTo(sm, 2, op2);
	CHK_EQUAL(0, tensorOpCheck(op1), errs);
	CHK_EQUAL(0, tensorOpCheck(op2), errs);


	elemOpDestroy(&sz);
	elemOpDestroy(&sp);
	elemOpDestroy(&sm);
	tensorOpDestroy(&op1);
	tensorOpDestroy(&op2);
	prodSpaceDestroy(&hTot);
	prodSpaceDestroy(&h1);
	return errs;
}

static int testTensorOpScale()
{
	int errs = 0;
	TensorOp a;
	ProdSpace h;
	double alpha = 2.0;

	h = prodSpaceCreate(0);
	tensorOpNull(h, &a);
	tensorOpScale(alpha, &a);
	CHK_EQUAL(a->numTerms, 0, errs);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	prodSpaceDestroy(&h);

	h = prodSpaceCreate(2);
	tensorOpNull(h, &a);
	tensorOpScale(alpha, &a);
	CHK_EQUAL(a->numTerms, 0, errs);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	prodSpaceDestroy(&h);

	h = prodSpaceCreate(2);
	tensorOpIdentity(h, &a);
	tensorOpScale(alpha, &a);
	CHK_EQUAL(a->numTerms, 1, errs);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	prodSpaceDestroy(&h);

	return errs;
}

static int testTensorOpCheck()
{
	int errs = 0;
	return errs;
}

static int testTensorOpKron()
{
	int errs = 0;
	TensorOp a;
	ProdSpace h;

	h = prodSpaceCreate(0);
	tensorOpNull(h, &a);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	tensorOpIdentity(h, &a);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	prodSpaceDestroy(&h);

	h = prodSpaceCreate(30);
	tensorOpNull(h, &a);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	tensorOpIdentity(h, &a);
	CHK_EQUAL(tensorOpCheck(a), 0, errs);
	tensorOpDestroy(&a);
	prodSpaceDestroy(&h);

	return errs;
}

int tensorOpTest()
{
	int errs = 0;
	errs += testTensorOpNull();
	errs += testTensorOpIdentity();
	errs += testTensorOpAddTo();
	errs += testTensorOpAddScaledTo();
	errs += testFindEmbedding();
	errs += testGatherIthEmbedding();
	errs += testTensorOpMul();
	errs += testTensorOpPlus();
	errs += testTensorOpScale();
	errs += testTensorOpCheck();
	errs += testTensorOpKron();
	return errs;
}
