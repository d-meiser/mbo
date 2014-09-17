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
		     struct MboAmplitude *x, struct MboAmplitude *y,
		     MboGlobInd rmin, MboGlobInd rmax)
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
				    x + entries[e].n * blockSizeAfter,
				    y + entries[e].m * blockSizeAfter, rmin,
				    rmax);
			}
			x += blockSizeAfter * (MboGlobInd)dims[nextI];
			y += blockSizeAfter * (MboGlobInd)dims[nextI];
		}
	} else {
		for (n = 0; n < blockSizeAfter; ++n) {
			y[n].re += alpha.re * x[n].re - alpha.im * x[n].im;
			y[n].im += alpha.re * x[n].im + alpha.im * x[n].re;
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
			   int numFactors, struct Embedding *embeddings,
			   MboGlobInd rmin, MboGlobInd rmax, int *I, int *J,
			   struct MboAmplitude *A, MboGlobInd offsetR,
			   MboGlobInd offsetC, int *numInserted)
{
	MboGlobInd blockSizeBefore, n, nStart, nEnd;
	struct MboAmplitude tmp;
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
				tmp.re = alpha.re * entries[e].val.re -
					 alpha.im * entries[e].val.im;
				tmp.im = alpha.re * entries[e].val.im +
					 alpha.im * entries[e].val.re;
				embeddingSparseMatrix(
				    nextI + 1, numSpaces, dims, blockSizeAfter,
				    tmp, numFactors - 1, embeddings + 1, rmin,
				    rmax,
				    I, J, A,
				    offsetR + entries[e].m * blockSizeAfter,
				    offsetC + entries[e].n * blockSizeAfter,
				    numInserted);
			}
			offsetR += blockSizeAfter * (MboGlobInd)dims[nextI];
			offsetC += blockSizeAfter * (MboGlobInd)dims[nextI];
		}
	} else {
		nStart = offsetR - rmin;
		if (nStart < 0) nStart = 0;
		nEnd = offsetR + blockSizeAfter - rmin;
		if (nEnd > rmax - rmin) nEnd = rmax - rmin;
		for (n = nStart; n < nEnd; ++n) {
			J[I[n] + numInserted[n]] = rmin + offsetC - offsetR + n;
			A[I[n] + numInserted[n]] = alpha;
			++numInserted[n];
		}
	}
}

void embeddingDiagonal(int i, int numSpaces, MboLocInd *dims,
		       MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
		       int numFactors, struct Embedding *embeddings,
		       MboGlobInd rmin, MboGlobInd rmax,
		       struct MboAmplitude *diag, MboGlobInd offsetR,
		       MboGlobInd offsetC)
{
	MboGlobInd blockSizeBefore, n, nStart, nEnd;
	struct MboAmplitude tmp;
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
				if (entries[e].m != entries[e].n) continue;
				tmp.re = alpha.re * entries[e].val.re -
					 alpha.im * entries[e].val.im;
				tmp.im = alpha.re * entries[e].val.im +
					 alpha.im * entries[e].val.re;
				embeddingDiagonal(
				    nextI + 1, numSpaces, dims, blockSizeAfter,
				    tmp, numFactors - 1, embeddings + 1, rmin,
				    rmax, diag,
				    offsetR + entries[e].m * blockSizeAfter,
				    offsetC + entries[e].n * blockSizeAfter);
			}
			offsetR += blockSizeAfter * (MboGlobInd)dims[nextI];
			offsetC += blockSizeAfter * (MboGlobInd)dims[nextI];
		}
	} else {
		nStart = offsetR - rmin;
		if (nStart < 0) nStart = 0;
		nEnd = offsetR + blockSizeAfter - rmin;
		if (nEnd > rmax - rmin) nEnd = rmax - rmin;
		for (n = nStart; n < nEnd; ++n) {
			diag[n] = alpha;
		}
	}
}

void embeddingDeleteDiagonal(struct Embedding *embedding)
{
	int e, numEntries;
	struct MboNonZeroEntry *entries;

	numEntries = mboElemOpNumEntries(embedding->op);
	e = 0;
	while (e < numEntries) {
		entries = mboElemOpGetEntries(embedding->op);
		if (entries[e].m == entries[e].n) {
			mboElemOpDeleteEntry(embedding->op, e);
			numEntries = mboElemOpNumEntries(embedding->op);
		} else {
			++e;
		}
	}
}
