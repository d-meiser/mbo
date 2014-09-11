#include <stdlib.h>
#include <string.h>

#include <MboTensorOp.h>
#include <MboIndices.h>
#include <MboAmplitude.h>
#include <MboNonZeroEntry.h>

#include <Embedding.h>
#include <SimpleTOp.h>
#include <Utilities.h>

/**
   Data structure for tensor product operators.
   */
struct MboTensorOp_t
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

MboProdSpace mboTensorOpGetSpace(MboTensorOp op)
{
	return op->space;
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

void mboTensorOpAddScaledTo(struct MboAmplitude *alpha, MboElemOp a, int i,
			    MboTensorOp op)
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

static void mboTensorOpKronTwo(MboTensorOp a, MboTensorOp b, MboTensorOp *c)
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

MBO_STATUS mboTensorOpKron(int n, MboTensorOp *ops, MboTensorOp *c)
{
	MboProdSpace h;
	int i, spacesEqual;
	MboTensorOp a, b;

	h = mboProdSpaceCreate(0);
	/* Below we're multiplying the spaces together starting at the right,
	 * but the operators are listed from left to right.  Hence we have to
	 * assemble here starting from the right to check the product space
	 * equality */
	for (i = n - 1; i >= 0; --i) {
		mboProdSpaceMul(ops[i]->space, &h);
	}
	spacesEqual = mboProdSpaceEqual(h, (*c)->space);
	mboProdSpaceDestroy(&h);
	if (!spacesEqual) return MBO_SPACE_MISMATCH;

	/* fold the operators together using mboTensorOpKronTwo for pairwise
	 * multiplication */
	h = mboProdSpaceCreate(0);
	mboTensorOpIdentity(h, &b);
	for (i = 0; i < n; ++i) {
		mboProdSpaceMul(ops[i]->space, &h);
		mboTensorOpNull(h, &a);
		mboTensorOpKronTwo(b, ops[i], &a);
		mboTensorOpDestroy(&b);
		b = a;
	}
	mboTensorOpPlus(b, c);
	mboTensorOpDestroy(&b);
	mboProdSpaceDestroy(&h);
	return MBO_SUCCESS;
}


MBO_STATUS mboTensorOpMatVec(struct MboAmplitude *alpha, MboTensorOp a,
			     MboVec x, struct MboAmplitude *beta, MboVec y)
{
	int i;
	MBO_STATUS err;
	MboVec tmp;
	struct MboAmplitude zero;

	if (mboProdSpaceDim(a->space) != mboVecDim(x) ||
	    mboProdSpaceDim(a->space) != mboVecDim(y)) {
		return MBO_DIMENSIONS_MISMATCH;
	}
	mboVecCreate(mboVecDim(y), &tmp);
	zero.re = 0;
	zero.im = 0;
	mboVecSet(&zero, tmp);
	for (i = 0; i < a->numTerms; ++i) {
		err = applySimpleTOp(a->space, alpha, a->sum + i, x, tmp);
		if (err != MBO_SUCCESS) return err;
	}
	mboVecAXPY(beta, y, tmp);
	mboVecSwap(tmp, y);
	mboVecDestroy(&tmp);
	return MBO_SUCCESS;
}

double mboTensorOpFlops(MboTensorOp a)
{
	int i, numSpaces;
	MboLocInd *dims;
	double flops = 0;

	dims = malloc(mboProdSpaceSize(a->space) * sizeof(*dims));
	numSpaces = mboProdSpaceSize(a->space);
	mboProdSpaceGetDims(a->space, numSpaces, dims);
	for (i = 0; i < a->numTerms; ++i) {
		flops += flopsSimpleTOp(numSpaces, dims, a->sum + i);
	}
	free(dims);
	return flops;
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

static void zeroArray(MboGlobInd n, struct MboAmplitude *array)
{
	MboGlobInd i;
	for (i = 0; i < n; ++i) {
		array[i].re = 0;
		array[i].im = 0;
	}
}

void mboTensorOpDenseMatrix(MboTensorOp a, struct MboAmplitude *mat)
{
	int i;
	MboGlobInd dim;
	dim = mboProdSpaceDim(a->space);

	zeroArray(dim * dim, mat);
	for (i = 0; i < a->numTerms; ++i) {
		simpleTOpDenseMatrix(a->space, a->sum + i, mat);
	}
}

void mboTensorOpRowOffsets(MboTensorOp a, MboGlobInd rmin, MboGlobInd rmax,
			   int *i)
{
	int j, r;
	for (r = 0; r < rmax - rmin; ++r) {
		i[r] = 0;
	}
	for (j = 0; j < a->numTerms; ++j) {
		simpleTOpGetNonZerosPerRow(a->space, a->sum + j, rmin, rmax,
					   i + 1);
	}
	i[0] = 0;
	for (j = 0; j < rmax - rmin; ++j) {
		i[j + 1] += i[j];
	}
}
