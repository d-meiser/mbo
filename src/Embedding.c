/*
Copyright 2014 Dominic Meiser

This file is part of mbo.

mbo is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

mbo is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with mbo.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <Embedding.h>
#include <stdlib.h>
#include <assert.h>
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

static MboGlobInd min(MboGlobInd a, MboGlobInd b)
{
	if (a < b) {
		return a;
	} else {
		return b;
	}
}

static MboGlobInd max(MboGlobInd a, MboGlobInd b)
{
	if (a > b) {
		return a;
	} else {
		return b;
	}
}

void applyLeafMask(struct MboAmplitude alpha, struct MboAmplitude *x,
		      struct MboAmplitude *y, const struct Tile *tile,
		      const struct Tile *mask)
{
	MboGlobInd rMin, cMin, rMax, r;
	struct Tile intersection = tileIntersection(*tile, *mask);

	rMin = tile->rmin + max(intersection.rmin - tile->rmin,
				intersection.cmin - tile->cmin);
	cMin = tile->cmin + max(intersection.rmin - tile->rmin,
				intersection.cmin - tile->cmin);
	rMax = tile->rmin + min(intersection.rmax - tile->rmin,
				intersection.cmax - tile->cmin);

	x -= rMin - max(cMin - mask->cmin, 0);
	y -= rMin - max(rMin - mask->rmin, 0);
	for (r = rMin; r < rMax; ++r) {
		y[r].re += alpha.re * x[r].re - alpha.im * x[r].im;
		y[r].im += alpha.re * x[r].im + alpha.im * x[r].re;
	}
}

static void applyLeaf0(struct MboAmplitude alpha, struct MboAmplitude *x,
		       struct MboAmplitude *y, MboGlobInd blockSize)
{
	MboGlobInd n;
	for (n = 0; n < blockSize; ++n) {
		y[n].re += alpha.re * x[n].re - alpha.im * x[n].im;
		y[n].im += alpha.re * x[n].im + alpha.im * x[n].re;
	}
}

static void applyLeaf1(struct MboAmplitude alpha, struct MboAmplitude *x,
		       struct MboAmplitude *y, MboGlobInd blockSizeBefore,
		       MboGlobInd blockSizeAfterAndSelf,
		       MboGlobInd blockSizeAfter)
{
	MboGlobInd m, n;
	for (m = 0; m < blockSizeBefore; ++m) {
		for (n = 0; n < blockSizeAfter; ++n) {
			y[n].re += alpha.re * x[n].re - alpha.im * x[n].im;
			y[n].im += alpha.re * x[n].im + alpha.im * x[n].re;
		}
		x += blockSizeAfterAndSelf;
		y += blockSizeAfterAndSelf;
	}
}

static void applyLeaf2(struct MboAmplitude alpha, struct MboAmplitude *x,
		       struct MboAmplitude *y, MboGlobInd M, MboGlobInd N,
		       MboGlobInd O, MboGlobInd P, MboGlobInd Q)
{
	MboGlobInd m, n, o;
	struct MboAmplitude *xtmp, *ytmp;
	for (m = 0; m < M; ++m) {
		xtmp = x;
		ytmp = y;
		for (n = 0; n < O; ++n) {
			for (o = 0; o < Q; ++o) {
				ytmp[o].re += alpha.re * xtmp[o].re -
					      alpha.im * xtmp[o].im;
				ytmp[o].im += alpha.re * xtmp[o].im +
					      alpha.im * xtmp[o].re;
			}
			xtmp += P;
			ytmp += P;
		}
		x += N;
		y += N;
	}
}

void applyEmbeddingsMask(int i, int numSpaces, MboLocInd *dims,
		     MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
		     int numFactors, struct Embedding *embeddings,
		     struct MboAmplitude *x, struct MboAmplitude *y,
		     struct Tile tile, const struct Tile *mask)
{
	int nextI, e;
	MboGlobInd blockSizeBefore, n, numTiles;
	struct MboNonZeroEntry *entries;
	struct MboAmplitude tmp;
	struct Tile myTile, childTile, quotient;

	if (numFactors > 0) {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);

		myTile.rmin = tile.rmin;
		myTile.rmax = myTile.rmin + blockSizeAfter * dims[nextI];
		myTile.cmin = tile.cmin;
		myTile.cmax = myTile.cmin + blockSizeAfter * dims[nextI];
		quotient = tileDivide(*mask, myTile);
		if (quotient.rmin > quotient.cmin) {
			tileAdvance(quotient.rmin, &myTile);
		} else {
			tileAdvance(quotient.cmin, &myTile);
		}
		numTiles = numTilesContained(*mask, myTile);
		for (e = 0; e < mboElemOpNumEntries(embeddings->op); ++e) {
			tmp.re = alpha.re * entries[e].val.re -
				 alpha.im * entries[e].val.im;
			tmp.im = alpha.re * entries[e].val.im +
				 alpha.im * entries[e].val.re;

			childTile.rmin =
			    myTile.rmin + entries[e].m * blockSizeAfter;
			childTile.rmax = childTile.rmin + blockSizeAfter;
			childTile.cmin =
			    myTile.cmin + entries[e].n * blockSizeAfter;
			childTile.cmax = childTile.cmin + blockSizeAfter;

			for (n = 0; n < numTiles; ++n) {
				applyEmbeddingsMask(nextI + 1, numSpaces, dims,
						blockSizeAfter, tmp,
						numFactors - 1, embeddings + 1,
						x, y, childTile, mask);
				childTile.rmin += myTile.rmax - myTile.rmin;
				childTile.rmax += myTile.rmax - myTile.rmin;
				childTile.cmin += myTile.cmax - myTile.cmin;
				childTile.cmax += myTile.cmax - myTile.cmin;
			}
		}
	} else {
		applyLeafMask(alpha, x, y, &tile, mask);
	}
}

void applyEmbeddings(int i, int numSpaces, MboLocInd *dims,
		     MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
		     int numFactors, struct Embedding *embeddings,
		     struct MboAmplitude *x, struct MboAmplitude *y)
{
	int nextI, e, ep;
	MboGlobInd blockSizeBefore, n, M, N, O, P, Q;
	struct MboNonZeroEntry *entries, *entriesB;
	struct MboAmplitude tmp, tmpp, *xtmp, *ytmp;

	if (numFactors == 0) {
		applyLeaf0(alpha, x, y, blockSizeAfter);
	} else if (numFactors == 1) {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);
		for (e = 0; e < mboElemOpNumEntries(embeddings->op); ++e) {
			tmp.re = alpha.re * entries[e].val.re -
				 alpha.im * entries[e].val.im;
			tmp.im = alpha.re * entries[e].val.im +
				 alpha.im * entries[e].val.re;
			xtmp = x + entries[e].n * blockSizeAfter;
			ytmp = y + entries[e].m * blockSizeAfter;
			applyLeaf1(tmp, xtmp, ytmp, blockSizeBefore,
				   blockSizeAfter * (MboGlobInd)dims[nextI],
				   blockSizeAfter);
		}
	} else if (numFactors == 2) {
		nextI = embeddings[0].i;
		M = computeBlockSize(nextI - i, dims + i);
		N = blockSizeAfter / M;
		O = computeBlockSize(embeddings[1].i - nextI - 1,
				     dims + nextI + 1);
		nextI = embeddings[1].i;
		P = computeBlockSize(numSpaces - nextI, dims + nextI);
		Q = computeBlockSize(numSpaces - nextI - 1, dims + nextI + 1);
		entries = mboElemOpGetEntries(embeddings[0].op);
		entriesB = mboElemOpGetEntries(embeddings[1].op);
		for (e = 0; e < mboElemOpNumEntries(embeddings[0].op); ++e) {
			tmp.re = alpha.re * entries[e].val.re -
				 alpha.im * entries[e].val.im;
			tmp.im = alpha.re * entries[e].val.im +
				 alpha.im * entries[e].val.re;
			for (ep = 0; ep < mboElemOpNumEntries(embeddings[1].op);
			     ++ep) {
				tmpp.re = tmp.re * entriesB[ep].val.re -
					  tmp.im * entriesB[ep].val.im;
				tmpp.im = tmp.re * entriesB[ep].val.im +
					  tmp.im * entriesB[ep].val.re;
				xtmp =
				    x + entries[e].n * O * P + entriesB[ep].n * Q;
				ytmp =
				    y + entries[e].m * O * P + entriesB[ep].m * Q;
				applyLeaf2(tmpp, xtmp, ytmp, M, N, O, P, Q);
			}
		}

	} else {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);

		for (e = 0; e < mboElemOpNumEntries(embeddings->op); ++e) {
			tmp.re = alpha.re * entries[e].val.re -
				 alpha.im * entries[e].val.im;
			tmp.im = alpha.re * entries[e].val.im +
				 alpha.im * entries[e].val.re;
			xtmp = x + entries[e].n * blockSizeAfter;
			ytmp = y + entries[e].m * blockSizeAfter;
			for (n = 0; n < blockSizeBefore; ++n) {
				applyEmbeddings(nextI + 1, numSpaces, dims,
						blockSizeAfter, tmp,
						numFactors - 1, embeddings + 1,
						xtmp, ytmp);
				xtmp += blockSizeAfter * (MboGlobInd)dims[nextI];
				ytmp += blockSizeAfter * (MboGlobInd)dims[nextI];
			}
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

MboGlobInd embeddingsDistanceFromDiagonal(int i, int numSpaces, MboLocInd *dims,
					  MboGlobInd blockSizeAfter,
					  int numFactors,
					  struct Embedding *embeddings)
{
	int nextI, e;
	MboGlobInd blockSizeBefore, meanDistance;
	struct MboNonZeroEntry *entries;

	if (numFactors > 0) {
		nextI = embeddings->i;
		blockSizeBefore = computeBlockSize(nextI - i, dims + i);
		blockSizeAfter /= (blockSizeBefore * (MboGlobInd)dims[nextI]);
		entries = mboElemOpGetEntries(embeddings->op);
		meanDistance = 0;
		for (e = 0; e < mboElemOpNumEntries(embeddings->op); ++e) {
			meanDistance +=
			    (entries[e].m - entries[e].n) * blockSizeAfter;
			meanDistance += embeddingsDistanceFromDiagonal(
			    nextI + 1, numSpaces, dims, blockSizeAfter,
			    numFactors - 1, embeddings + 1);
		}
		meanDistance /= mboElemOpNumEntries(embeddings->op);
		return meanDistance;
	} else {
		return 0;
	}
}
