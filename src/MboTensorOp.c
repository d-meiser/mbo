#include <stdlib.h>
#include <string.h>

#include <MboTensorOp.h>
#include <MboIndices.h>
#include <MboAmplitude.h>
#include <MboNonZeroEntry.h>

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
/**
 * @brief Merge embeddings for each slot
 * Effectively calls gatherIthEmbedding for each index in the *embeddings array.
 * */
void gatherAllEmbeddings(int *numEmbeddings, struct Embedding **embeddings);
/**
 * @brief sort embeddings in ascending order according to i */
static void sortEmbeddings(int numEmbeddings, struct Embedding *a);
static void destroySimpleTOp(struct SimpleTOp *term);
static void multiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);
static void kronSimpleTOps(struct SimpleTOp *a, int numSpacesInA,
			   struct SimpleTOp *b, struct SimpleTOp *c);
static void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);
static void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h,
			   struct SimpleTOp *op);
static int checkSimpleTOp(struct SimpleTOp *sa);
static MBO_STATUS applySimpleTOp(MboProdSpace h, struct MboAmplitude *alpha,
				 struct SimpleTOp *a, MboVec x, MboVec y);
static MboGlobInd computeBlockSize(int N, MboLocInd *dims);
static void applyEmbeddings(int i, int numSpaces, MboLocInd *dims,
			    MboGlobInd blockSize, struct MboAmplitude alpha,
			    int numFactors, struct Embedding *embeddings,
			    struct MboAmplitude *xarr,
			    struct MboAmplitude *yarr);
static double flopsSimpleTOp(int numSpaces, MboLocInd *dims,
			     struct SimpleTOp *op);

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
				mboElemOpCreate(
				    &b->embeddings[b->numFactors].op);
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

void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h,
		    struct SimpleTOp *op)
{
	MboLocInd d;
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

MBO_STATUS applySimpleTOp(MboProdSpace h, struct MboAmplitude *alpha,
			  struct SimpleTOp *a, MboVec x, MboVec y)
{
	int numSpaces;
	MboLocInd *dims;
	MboGlobInd blockSize;
	struct MboAmplitude *xarr, *yarr;

	dims = malloc(mboProdSpaceSize(h) * sizeof(*dims));
	gatherAllEmbeddings(&a->numFactors, &a->embeddings);
	sortEmbeddings(a->numFactors, a->embeddings);

	numSpaces = mboProdSpaceSize(h);
	blockSize = mboProdSpaceDim(h);
	mboProdSpaceGetDims(h, mboProdSpaceSize(h), dims);

	mboVecGetViewR(x, &xarr);
	mboVecGetViewRW(y, &yarr);
	applyEmbeddings(0, numSpaces, dims, blockSize, *alpha, a->numFactors,
			a->embeddings, xarr, yarr);
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);

	free(dims);

	return MBO_SUCCESS;
}

void applyEmbeddings(int i, int numSpaces, MboLocInd *dims,
		     MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
		     int numFactors, struct Embedding *embeddings,
		     struct MboAmplitude *xarr, struct MboAmplitude *yarr)
{
	int nextI, e;
	MboGlobInd blockSizeBefore, n;
	struct MboNonZeroEntry *entries;
	struct MboAmplitude tmp;

	if (numFactors > 0) {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);
		for (n = 0; n < blockSizeBefore; ++n) {
			for (e = 0; e < mboElemOpNumEntries(embeddings->op);
					++e) {
				tmp.re = alpha.re * entries[e].val.re -
					       alpha.im * entries[e].val.im;
				tmp.im = alpha.re * entries[e].val.im +
					       alpha.im * entries[e].val.re;
				applyEmbeddings(
				    nextI + 1, numSpaces, dims, blockSizeAfter,
				    tmp, numFactors - 1, embeddings + 1,
				    xarr + entries[e].n * blockSizeAfter,
				    yarr + entries[e].m * blockSizeAfter);
			}
			xarr += blockSizeAfter * (MboGlobInd)dims[nextI];
			yarr += blockSizeAfter * (MboGlobInd)dims[nextI];
		}
	} else {
		for (n = 0; n < blockSizeAfter; ++n) {
			yarr[n].re +=
			    alpha.re * xarr[n].re - alpha.im * xarr[n].im;
			yarr[n].im +=
			    alpha.re * xarr[n].im + alpha.im * xarr[n].re;
		}
	}
}

void gatherAllEmbeddings(int *numEmbeddings, struct Embedding **embeddings)
{
	int *is, nInitial = *numEmbeddings, i;
	is = malloc(*numEmbeddings * sizeof(*is));
	for (i = 0; i < nInitial; ++i) {
		is[i] = (*embeddings)[i].i;
	}
	for (i = 0; i < nInitial; ++i) {
		gatherIthEmbedding(is[i], numEmbeddings, embeddings);
	}
	free(is);
}

static int embeddingCmp(const void *p1, const void *p2)
{
	struct Embedding *e1 = (struct Embedding *)p1;
	struct Embedding *e2 = (struct Embedding *)p2;
	return e1->i - e2->i;
}

void sortEmbeddings(int numEmbeddings, struct Embedding *embeddings)
{
	qsort(embeddings, numEmbeddings, sizeof(*embeddings), embeddingCmp);
}

MboGlobInd computeBlockSize(int N, MboLocInd *dims)
{
	int i;
	MboGlobInd blockSize = 1;
	for (i = 0; i < N; ++i) {
		blockSize *= (MboGlobInd)dims[i];
	}
	return blockSize;
}

double flopsSimpleTOp(int numSpaces, MboLocInd *dims, struct SimpleTOp *a)
{
	int i, j;
	double flops;

	gatherAllEmbeddings(&a->numFactors, &a->embeddings);
	sortEmbeddings(a->numFactors, a->embeddings);

	flops = 1;
	for (i = 0; i < numSpaces; ++i) {
		j = gatherIthEmbedding(i, &a->numFactors, &a->embeddings);
		if (j < 0) {
			flops *= dims[i];
		} else {
			flops *= mboElemOpNumEntries(a->embeddings[j].op);
		}
	}
	return flops * 8.0;
}

