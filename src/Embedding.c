#include <Embedding.h>
#include <stdlib.h>
#include <Utilities.h>
#include <MboNonZeroEntry.h>

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

void embeddingDenseMatrix(int i, int numSpaces, MboLocInd *dims,
			  MboGlobInd blockSizeAfter, MboGlobInd dim,
			  struct MboAmplitude alpha, int numFactors,
			  struct Embedding *embeddings,
			  struct MboAmplitude *mat)
{
	MboGlobInd blockSizeBefore, n;
	struct MboNonZeroEntry *entries;
	int nextI, e;
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
				embeddingDenseMatrix(
				    nextI + 1, numSpaces, dims, blockSizeAfter,
				    dim, tmp, numFactors - 1, embeddings + 1,
				    mat + entries[e].n * blockSizeAfter +
					entries[e].m * blockSizeAfter * dim);
			}
			mat += blockSizeAfter * (MboGlobInd)dims[nextI] *
			       (dim + 1);
		}
	} else {
		for (n = 0; n < blockSizeAfter; ++n) {
			mat[0].re += alpha.re;
			mat[0].im += alpha.im;
			mat += dim + 1;
		}
	}
}

void embeddingNonZeros(int i, int numSpaces, MboLocInd *dims,
		       MboGlobInd blockSizeAfter, int numFactors,
		       struct Embedding *embeddings, MboGlobInd rmin,
		       MboGlobInd rmax, MboGlobInd offset, int *nnz)
{
	MboGlobInd blockSizeBefore, n, nStart, nEnd;
	struct MboNonZeroEntry *entries;
	int nextI, e;

	if (numFactors > 0) {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);
		for (n = 0; n < blockSizeBefore; ++n) {
			for (e = 0; e < mboElemOpNumEntries(embeddings->op);
			     ++e) {
				embeddingNonZeros(
				    nextI + 1, numSpaces, dims, blockSizeAfter,
				    numFactors - 1, embeddings + 1, rmin, rmax,
				    offset + entries[e].m * blockSizeAfter,
				    nnz);
			}
			offset += blockSizeAfter * (MboGlobInd)dims[nextI];
		}
	} else {
		nStart = offset - rmin;
		if (nStart < 0) nStart = 0;
		nEnd = offset + blockSizeAfter - rmin;
		if (nEnd > rmax - rmin) nEnd = rmax - rmin;
		for (n = nStart; n < nEnd; ++n) {
			++nnz[n];
		}
	}
}

void embeddingSparseMatrix(int i, int numSpaces, MboLocInd *dims,
			   MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
			   int numFactors, struct Embedding *embeddings, int *I,
			   int *J, struct MboAmplitude *A, int *numInserted)
{
}
