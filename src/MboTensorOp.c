#include <stdlib.h>

#include <MboTensorOp.h>
#include <MboAmplitude.h>

/**
 * @brief Elementary embedding of an operator into a product space.
 * */
struct Embedding
{
	int i;
	MboElemOp op;
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
static void kronSimpleTOps(struct SimpleTOp *a, int numSpacesInA,
			   struct SimpleTOp *b, struct SimpleTOp *c);
static void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);
static void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h,
			   struct SimpleTOp *op);
static int checkSimpleTOp(struct SimpleTOp *sa);

/**
   Data structure for tensor product operators.
   */
struct MboTensorOp
{
	/* Product space on which the operator is defined. */
	MboProdSpace space;
	/* Number of terms (SimpleTOp's) making up the operator */
	int numTerms;
	/* Array of terms in sum */
	struct SimpleTOp *sum;
};

void mboTensorOpNull(MboProdSpace h, MboTensorOp *op)
{
	MboTensorOp a = malloc(sizeof(*a));
	a->space = mboProdSpaceCopy(h);
	a->numTerms = 0;
	a->sum = 0;
	*op = a;
}

void mboTensorOpIdentity(MboProdSpace h, MboTensorOp *op)
{
	MboTensorOp a = malloc(sizeof(*a));
	a->space = mboProdSpaceCopy(h);
	a->numTerms = 1;
	a->sum = malloc(sizeof(*a->sum));
	a->sum[0].numFactors = 0;
	a->sum[0].embeddings = 0;
	*op = a;
}

void mboTensorOpDestroy(MboTensorOp *op)
{
	int i;
	for (i = 0; i < (*op)->numTerms; ++i) {
		destroySimpleTOp(&(*op)->sum[i]);
	}
	free((*op)->sum);
	mboProdSpaceDestroy(&(*op)->space);
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

void mboTensorOpAddTo(MboElemOp a, int i, MboTensorOp op)
{
	struct SimpleTOp *newTerm;
	struct Embedding *aEmbedded;

	op->sum = realloc(op->sum, (op->numTerms + 1) * sizeof(*op->sum));
	newTerm = op->sum + op->numTerms;
	newTerm->numFactors = 1;
	newTerm->embeddings = malloc(sizeof(struct Embedding));
	aEmbedded = newTerm->embeddings;
	aEmbedded->i = i;
	aEmbedded->op = mboElemOpCopy(a);
	++op->numTerms;
}

void mboTensorOpAddScaledTo(struct MboAmplitude *alpha, MboElemOp a, int i, MboTensorOp op)
{
	struct SimpleTOp *newTerm;
	struct Embedding *aEmbedded;

	op->sum = realloc(op->sum, (op->numTerms + 1) * sizeof(*op->sum));
	newTerm = op->sum + op->numTerms;
	newTerm->numFactors = 1;
	newTerm->embeddings = malloc(sizeof(struct Embedding));
	aEmbedded = newTerm->embeddings;
	aEmbedded->i = i;
	aEmbedded->op = mboElemOpCopy(a);
	mboElemOpScale(alpha, aEmbedded->op);
	++op->numTerms;
}

void mboTensorOpMul(MboTensorOp a, MboTensorOp b, MboTensorOp *c)
{
	int i, j;
	int N;
	int numNewTerms = a->numTerms * b->numTerms;
	struct SimpleTOp *cFillPtr;

	if (numNewTerms) {
		(*c)->sum = realloc((*c)->sum, ((*c)->numTerms + numNewTerms) *
						   (sizeof(*(*c)->sum)));
	}

	N = mboProdSpaceSize(b->space);
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

void mboTensorOpPlus(MboTensorOp a, MboTensorOp *b)
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

void mboTensorOpScale(struct MboAmplitude *alpha, MboTensorOp *a)
{
	int i;
	for (i = 0; i < (*a)->numTerms; ++i) {
		scaleSimpleTOp(alpha, (*a)->space, (*a)->sum + i);
	}
}

void mboTensorOpKron(MboTensorOp a, MboTensorOp b, MboTensorOp *c)
{
	int i, j, numNewTerms = a->numTerms * b->numTerms;
	struct SimpleTOp *cFillPtr;
	(*c)->sum = realloc((*c)->sum, ((*c)->numTerms + numNewTerms) *
					   sizeof(*(*c)->sum));
	cFillPtr = (*c)->sum + (*c)->numTerms;
	for (i = 0; i < a->numTerms; ++i) {
		for (j = 0; j < b->numTerms; ++j) {
			cFillPtr->numFactors = 0;
			cFillPtr->embeddings = 0;
			kronSimpleTOps(a->sum + i, mboProdSpaceSize(a->space),
				       b->sum + j, cFillPtr);
			++cFillPtr;
		}
	}
	(*c)->numTerms += numNewTerms;
}

void kronSimpleTOps(struct SimpleTOp *a, int numSpacesInA, struct SimpleTOp *b,
		    struct SimpleTOp *c)
{
	int i;
	struct Embedding *cFillPtr;
	c->embeddings = realloc(
	    c->embeddings, (c->numFactors + a->numFactors + b->numFactors) *
			       sizeof(*c->embeddings));
	cFillPtr = c->embeddings + c->numFactors;
	for (i = 0; i < a->numFactors + b->numFactors; ++i) {
		mboElemOpCreate(&cFillPtr[i].op);
	}
	for (i = 0; i < a->numFactors; ++i) {
		copyEmbedding(cFillPtr + i, a->embeddings + i);
	}
	for (i = 0; i < b->numFactors; ++i) {
		copyEmbedding(cFillPtr + i + a->numFactors, b->embeddings + i);
		cFillPtr[a->numFactors + i].i += numSpacesInA;
	}
	c->numFactors += a->numFactors + b->numFactors;
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
				mboElemOpCreate(&b->embeddings[b->numFactors].op);
				copyEmbedding(b->embeddings + b->numFactors,
						a->embeddings + j);
				++b->numFactors;
			} else {
				/* Non-trivial operators in both slots; just
				   multiply them together. */
				mboElemOpMul(a->embeddings[j].op,
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
		mboElemOpCreate(&dest->embeddings[i].op);
		copyEmbedding(dest->embeddings + i, src->embeddings + i);
	}
}

void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h, struct SimpleTOp *op)
{
	int d;
	if (mboProdSpaceSize(h) == 0) return;
	if (op->numFactors == 0) {
		op->embeddings = malloc(sizeof(*op->embeddings));
		op->embeddings[0].i = 0;
		mboProdSpaceGetDims(h, 1, &d);
		op->embeddings[0].op = mboEye(d);
		op->numFactors += 1;
	}
	mboElemOpScale(alpha, op->embeddings[0].op);
}

void copyEmbedding(struct Embedding *dest, struct Embedding *src)
{
	dest->i = src->i;
	if (dest->op) {
		mboElemOpDestroy(&dest->op);
	}
	dest->op = mboElemOpCopy(src->op);
}

void destroyEmbedding(struct Embedding *e)
{
	mboElemOpDestroy(&e->op);
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
			mboElemOpPlus((*embeddings)[next].op,
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

int mboTensorOpCheck(MboTensorOp op)
{
	int i, errs = 0;
	if (op->numTerms < 0) ++errs;
	errs += mboProdSpaceCheck(op->space);
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
		errs += mboElemOpCheck(sa->embeddings[i].op);
	}
	return errs;
}

/*
 * Tests
 * */
#include "TestUtils.h"

#define EPS 1.0e-12

static int testMboTensorOpNull()
{
	int errs = 0;
	MboTensorOp op;
	MboProdSpace h = mboProdSpaceCreate(1);
	mboTensorOpNull(h, &op);
	CHK_TRUE(mboProdSpaceEqual(op->space, h), errs);
	CHK_EQUAL(op->numTerms, 0, errs);
	mboTensorOpDestroy(&op);
	mboProdSpaceDestroy(&h);
	return errs;
}

static int testMboTensorOpIdentity()
{
	int errs = 0;
	MboTensorOp op;
	MboProdSpace h = mboProdSpaceCreate(1);
	mboTensorOpIdentity(h, &op);
	CHK_TRUE(mboProdSpaceEqual(op->space, h), errs);
	CHK_EQUAL(op->numTerms, 1, errs);
	mboTensorOpDestroy(&op);
	mboProdSpaceDestroy(&h);
	return errs;
}

static int testMboTensorOpAddTo()
{
	int errs = 0;
	MboTensorOp op;
	MboElemOp eop;
	MboProdSpace h;
	struct MboAmplitude alpha;

	h = mboProdSpaceCreate(20);

	mboElemOpCreate(&eop);
	alpha.re = 5.4;
	alpha.im = 2.9;
	mboElemOpAddTo(14, 15, &alpha, &eop);

	mboTensorOpNull(h, &op);
	mboTensorOpAddTo(eop, 0, op);
	CHK_EQUAL(op->numTerms, 1, errs);
	CHK_TRUE(mboProdSpaceEqual(op->space, h), errs);
	CHK_TRUE(op->sum != 0, errs);
	CHK_EQUAL(op->sum[0].numFactors, 1, errs);
	CHK_TRUE(op->sum[0].embeddings != 0, errs);
	CHK_EQUAL(op->sum[0].embeddings[0].i, 0, errs);

	mboTensorOpDestroy(&op);
	mboElemOpDestroy(&eop);
	mboProdSpaceDestroy(&h);
	return errs;
}

static int testMboTensorOpAddScaledTo()
{
	int errs = 0;
	MboTensorOp op;
	MboElemOp eop;
	MboProdSpace h;
	struct MboAmplitude alpha, beta;

	alpha.re = 20;
	alpha.im = 30;
	beta.re = 15;
	beta.im = -19;

	h = mboProdSpaceCreate(20);

	mboElemOpCreate(&eop);
	mboElemOpAddTo(14, 15, &alpha, &eop);

	mboTensorOpNull(h, &op);
	mboTensorOpAddScaledTo(&alpha, eop, 0, op);
	CHK_EQUAL(op->numTerms, 1, errs);
	CHK_TRUE(mboProdSpaceEqual(op->space, h), errs);
	CHK_TRUE(op->sum != 0, errs);
	CHK_EQUAL(op->sum[0].numFactors, 1, errs);
	CHK_TRUE(op->sum[0].embeddings != 0, errs);
	CHK_EQUAL(op->sum[0].embeddings[0].i, 0, errs);

	mboTensorOpAddScaledTo(&beta, eop, 1, op);
	CHK_EQUAL(op->numTerms, 2, errs);

	mboTensorOpDestroy(&op);
	mboElemOpDestroy(&eop);
	mboProdSpaceDestroy(&h);
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
	mboElemOpCreate(&embeddings[0].op);
	embeddings[1].i = 4;
	mboElemOpCreate(&embeddings[1].op);
	embeddings[2].i = 0;
	mboElemOpCreate(&embeddings[2].op);
	embeddings[3].i = 3;
	mboElemOpCreate(&embeddings[3].op);
	embeddings[4].i = 0;
	mboElemOpCreate(&embeddings[4].op);

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
		mboElemOpDestroy(&embeddings[i].op);
	}

	free(embeddings);

	return errs;
}

static int testMboTensorOpMul()
{
	int errs = 0;
	int i;
	int N = 10;
	MboTensorOp op1, op2, op3;
	MboElemOp eop1, eop2;
	MboProdSpace h1, h2;
	struct MboAmplitude alpha, beta;

	alpha.re = 20;
	alpha.im = 30;
	beta.re = 15;
	beta.im = -19;

	h1 = mboProdSpaceCreate(5);
	h2 = mboProdSpaceCreate(0);
	for (i = 0; i < N; ++i) {
		mboProdSpaceMul(h1, &h2);
	}

	mboElemOpCreate(&eop1);
	mboElemOpAddTo(14, 15, &alpha, &eop1);
	mboElemOpAddTo(1, 5, &beta, &eop1);
	mboElemOpAddTo(2, 3, &beta, &eop1);

	mboElemOpCreate(&eop2);
	mboElemOpAddTo(8, 3, &beta, &eop2);
	mboElemOpAddTo(3, 4, &alpha, &eop2);

	mboTensorOpNull(h2, &op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboTensorOpNull(h2, &op1);
	mboTensorOpAddTo(eop1, 0, op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboTensorOpNull(h2, &op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpAddTo(eop1, 0, op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 0, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboTensorOpNull(h2, &op1);
	mboTensorOpAddTo(eop1, 0, op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpAddTo(eop2, 0, op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 1, errs);
	CHK_EQUAL(op3->sum[0].numFactors, 1, errs);
	CHK_EQUAL(op3->sum[0].embeddings[0].i, 0, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboTensorOpNull(h2, &op1);
	mboTensorOpAddTo(eop1, 0, op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpAddTo(eop2, 1, op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 1, errs);
	CHK_EQUAL(op3->sum[0].numFactors, 2, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboTensorOpNull(h2, &op1);
	mboTensorOpAddTo(eop1, 0, op1);
	mboTensorOpAddTo(eop1, 1, op1);
	mboTensorOpNull(h2, &op2);
	mboTensorOpAddTo(eop2, 1, op2);
	mboTensorOpAddTo(eop2, 2, op2);
	mboTensorOpNull(h2, &op3);
	mboTensorOpMul(op1, op2, &op3);
	CHK_EQUAL(op3->numTerms, 4, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboTensorOpDestroy(&op3);

	mboElemOpDestroy(&eop1);
	mboElemOpDestroy(&eop2);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	return errs;
}

static int testMboTensorOpPlus()
{
	int i, errs = 0, N = 5;
	MboProdSpace hTot, h1;
	MboTensorOp op1, op2;
	MboElemOp sz, sp, sm;

	hTot = mboProdSpaceCreate(0);
	h1 = mboProdSpaceCreate(2);
	for (i = 0; i < N; ++i) {
		mboProdSpaceMul(h1, &hTot);
	}

	sz = mboSigmaZ();
	sp = mboSigmaPlus();
	sm = mboSigmaMinus();

	mboTensorOpNull(hTot, &op1);
	mboTensorOpNull(hTot, &op2);
	mboTensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 0, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);

	mboTensorOpNull(hTot, &op1);
	mboTensorOpNull(hTot, &op2);
	mboTensorOpAddTo(sz, 1, op2);
	mboTensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 1, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);

	mboTensorOpNull(hTot, &op1);
	mboTensorOpAddTo(sz, 1, op1);
	mboTensorOpNull(hTot, &op2);
	mboTensorOpPlus(op1, &op2);
	CHK_EQUAL(op2->numTerms, 1, errs);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);

	mboTensorOpNull(hTot, &op1);
	mboTensorOpNull(hTot, &op2);
	mboTensorOpAddTo(sz, 0, op1);
	mboTensorOpAddTo(sz, 1, op1);
	mboTensorOpAddTo(sp, 0, op2);
	mboTensorOpAddTo(sp, 2, op2);
	mboTensorOpAddTo(sp, 3, op2);
	mboTensorOpAddTo(sm, 3, op2);
	mboTensorOpAddTo(sm, 2, op2);
	CHK_EQUAL(0, mboTensorOpCheck(op1), errs);
	CHK_EQUAL(0, mboTensorOpCheck(op2), errs);


	mboElemOpDestroy(&sz);
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sm);
	mboTensorOpDestroy(&op1);
	mboTensorOpDestroy(&op2);
	mboProdSpaceDestroy(&hTot);
	mboProdSpaceDestroy(&h1);
	return errs;
}

static int testMboTensorOpScale()
{
	int errs = 0;
	MboTensorOp a;
	MboProdSpace h;
	struct MboAmplitude alpha;

	alpha.re = 2.0;
	alpha.im = -55.55;

	h = mboProdSpaceCreate(0);
	mboTensorOpIdentity(h, &a);
	mboTensorOpScale(&alpha, &a);
	CHK_EQUAL(a->numTerms, 1, errs);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	h = mboProdSpaceCreate(2);
	mboTensorOpNull(h, &a);
	mboTensorOpScale(&alpha, &a);
	CHK_EQUAL(a->numTerms, 0, errs);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	h = mboProdSpaceCreate(2);
	mboTensorOpIdentity(h, &a);
	mboTensorOpScale(&alpha, &a);
	CHK_EQUAL(a->numTerms, 1, errs);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	return errs;
}

static int testMboTensorOpCheck()
{
	int errs = 0;
	MboTensorOp a;
	MboProdSpace h;

	h = mboProdSpaceCreate(0);
	mboTensorOpNull(h, &a);
	a->numTerms = -1;
	CHK_TRUE(mboTensorOpCheck(a) != 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	h = mboProdSpaceCreate(0);
	mboTensorOpNull(h, &a);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpIdentity(h, &a);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	h = mboProdSpaceCreate(30);
	mboTensorOpNull(h, &a);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpIdentity(h, &a);
	CHK_EQUAL(mboTensorOpCheck(a), 0, errs);
	mboTensorOpDestroy(&a);
	mboProdSpaceDestroy(&h);

	return errs;
}

static int testMboTensorOpKron()
{
	int errs = 0;
	MboProdSpace h1, h2, h3;
	MboTensorOp a, b, c;
	MboElemOp sz;
	MboElemOp sp;

	h1 = mboProdSpaceCreate(0);
	h2 = mboProdSpaceCreate(0);
	h3 = mboProdSpaceCreate(0);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h2, &h3);
	mboTensorOpNull(h1, &a);
	mboTensorOpNull(h2, &b);
	mboTensorOpNull(h3, &c);
	mboTensorOpKron(a, b, &c);
	CHK_EQUAL(mboTensorOpCheck(c), 0, errs);
	CHK_EQUAL(c->numTerms, 0, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpDestroy(&b);
	mboTensorOpDestroy(&c);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	h1 = mboProdSpaceCreate(0);
	h2 = mboProdSpaceCreate(0);
	h3 = mboProdSpaceCreate(0);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h2, &h3);
	mboTensorOpIdentity(h1, &a);
	mboTensorOpIdentity(h2, &b);
	mboTensorOpNull(h3, &c);
	mboTensorOpKron(a, b, &c);
	CHK_EQUAL(mboTensorOpCheck(c), 0, errs);
	CHK_EQUAL(c->numTerms, 1, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpDestroy(&b);
	mboTensorOpDestroy(&c);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	h1 = mboProdSpaceCreate(0);
	h2 = mboProdSpaceCreate(0);
	h3 = mboProdSpaceCreate(0);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h2, &h3);
	mboTensorOpIdentity(h1, &a);
	mboTensorOpIdentity(h2, &b);
	mboTensorOpIdentity(h3, &c);
	mboTensorOpKron(a, b, &c);
	CHK_EQUAL(mboTensorOpCheck(c), 0, errs);
	CHK_EQUAL(c->numTerms, 2, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpDestroy(&b);
	mboTensorOpDestroy(&c);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCreate(2);
	h3 = mboProdSpaceCreate(0);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h2, &h3);
	mboTensorOpIdentity(h1, &a);
	mboTensorOpIdentity(h2, &b);
	mboTensorOpIdentity(h3, &c);
	mboTensorOpKron(a, b, &c);
	CHK_EQUAL(mboTensorOpCheck(c), 0, errs);
	CHK_EQUAL(c->numTerms, 2, errs);
	mboTensorOpDestroy(&a);
	mboTensorOpDestroy(&b);
	mboTensorOpDestroy(&c);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCreate(2);
	h3 = mboProdSpaceCreate(0);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h2, &h3);
	mboTensorOpIdentity(h1, &a);
	sz = mboSigmaZ();
	mboTensorOpAddTo(sz, 0, a);
	mboTensorOpIdentity(h2, &b);
	sp = mboSigmaPlus();
	mboTensorOpAddTo(sp, 0, b);
	mboTensorOpIdentity(h3, &c);
	mboTensorOpKron(a, b, &c);
	CHK_EQUAL(mboTensorOpCheck(c), 0, errs);
	CHK_EQUAL(c->numTerms, 5, errs);
	mboElemOpDestroy(&sz);
	mboElemOpDestroy(&sp);
	mboTensorOpDestroy(&a);
	mboTensorOpDestroy(&b);
	mboTensorOpDestroy(&c);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	return errs;
}

static int testKronSimpleTOps()
{
	int errs = 0;
	struct SimpleTOp a, b, c;

	a.numFactors = 0;
	a.embeddings = 0;
	b.numFactors = 0;
	b.embeddings = 0;
	c.numFactors = 0;
	c.embeddings = 0;
	kronSimpleTOps(&a, 0, &b, &c);
	CHK_EQUAL(c.numFactors, 0, errs);
	destroySimpleTOp(&a);
	destroySimpleTOp(&b);
	destroySimpleTOp(&c);

	a.numFactors = 0;
	a.embeddings = 0;
	b.numFactors = 0;
	b.embeddings = 0;
	c.numFactors = 1;
	c.embeddings = malloc(sizeof(*c.embeddings));
	c.embeddings[0].i = 0;
	mboElemOpCreate(&c.embeddings[0].op);
	kronSimpleTOps(&a, 0, &b, &c);
	CHK_EQUAL(c.numFactors, 1, errs);
	destroySimpleTOp(&a);
	destroySimpleTOp(&b);
	destroySimpleTOp(&c);

	a.numFactors = 1;
	a.embeddings = malloc(sizeof(*a.embeddings));
	a.embeddings[0].i = 0;
	mboElemOpCreate(&a.embeddings[0].op);
	b.numFactors = 0;
	b.embeddings = 0;
	c.numFactors = 0;
	c.embeddings = 0;
	kronSimpleTOps(&a, 1, &b, &c);
	CHK_EQUAL(c.numFactors, 1, errs);
	destroySimpleTOp(&a);
	destroySimpleTOp(&b);
	destroySimpleTOp(&c);

	a.numFactors = 1;
	a.embeddings = malloc(sizeof(*a.embeddings));
	a.embeddings[0].i = 0;
	mboElemOpCreate(&a.embeddings[0].op);
	b.numFactors = 1;
	b.embeddings = malloc(sizeof(*b.embeddings));;
	b.embeddings[0].i = 0;
	mboElemOpCreate(&b.embeddings[0].op);
	c.numFactors = 0;
	c.embeddings = 0;
	kronSimpleTOps(&a, 1, &b, &c);
	CHK_EQUAL(c.numFactors, 2, errs);
	CHK_EQUAL(c.embeddings[0].i, 0, errs);
	CHK_EQUAL(c.embeddings[1].i, 1, errs);
	destroySimpleTOp(&a);
	destroySimpleTOp(&b);
	destroySimpleTOp(&c);
	return errs;
}

int mboTensorOpTest()
{
	int errs = 0;
	errs += testMboTensorOpNull();
	errs += testMboTensorOpIdentity();
	errs += testMboTensorOpAddTo();
	errs += testMboTensorOpAddScaledTo();
	errs += testFindEmbedding();
	errs += testGatherIthEmbedding();
	errs += testMboTensorOpMul();
	errs += testMboTensorOpPlus();
	errs += testMboTensorOpScale();
	errs += testMboTensorOpCheck();
	errs += testMboTensorOpKron();
	errs += testKronSimpleTOps();
	return errs;
}