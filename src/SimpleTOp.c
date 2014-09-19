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
#include <SimpleTOp.h>
#include <stdlib.h>

void destroySimpleTOp(struct SimpleTOp *term)
{
	int i;
	for (i = 0; i < term->numFactors; ++i) {
		destroyEmbedding(term->embeddings + i);
	}
	free(term->embeddings);
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

MBO_STATUS applySimpleTOp(MboProdSpace h, struct MboAmplitude alpha,
			  struct SimpleTOp *a, struct MboAmplitude *x,
			  struct MboAmplitude *y, MboGlobInd rmin,
			  MboGlobInd rmax)
{
	int numSpaces;
	MboLocInd *dims;
	MboGlobInd blockSize;

	dims = malloc(mboProdSpaceSize(h) * sizeof(*dims));
	gatherAllEmbeddings(&a->numFactors, &a->embeddings);
	sortEmbeddings(a->numFactors, a->embeddings);

	numSpaces = mboProdSpaceSize(h);
	blockSize = mboProdSpaceDim(h);
	mboProdSpaceGetDims(h, mboProdSpaceSize(h), dims);

	applyEmbeddings(0, numSpaces, dims, blockSize, alpha, a->numFactors,
			a->embeddings, x, y, rmin, rmax);

	free(dims);

	return MBO_SUCCESS;
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

void simpleTOpDenseMatrix(MboProdSpace h, struct SimpleTOp *simpleOp,
			  struct MboAmplitude *mat)
{
	MboGlobInd blockSize, dim;
	MboLocInd *dims;
	int numSpaces;
	struct MboAmplitude alpha;
	
	gatherAllEmbeddings(&simpleOp->numFactors, &simpleOp->embeddings);
	sortEmbeddings(simpleOp->numFactors, simpleOp->embeddings);

	blockSize = mboProdSpaceDim(h);
	dim = blockSize;
	numSpaces = mboProdSpaceSize(h);
	dims = malloc(numSpaces * sizeof(*dims));
	mboProdSpaceGetDims(h, numSpaces, dims);

	alpha.re = 1;
	alpha.im = 0;
	embeddingDenseMatrix(0, numSpaces, dims, blockSize, dim, alpha,
			simpleOp->numFactors, simpleOp->embeddings, mat);

	free(dims);
}

void simpleTOpGetNonZerosPerRow(MboProdSpace h, struct SimpleTOp *simpleOp,
				MboGlobInd rmin, MboGlobInd rmax, int *nnz)
{
	MboGlobInd blockSize, offset;
	MboLocInd *dims;
	int numSpaces;

	gatherAllEmbeddings(&simpleOp->numFactors, &simpleOp->embeddings);
	sortEmbeddings(simpleOp->numFactors, simpleOp->embeddings);

	blockSize = mboProdSpaceDim(h);
	numSpaces = mboProdSpaceSize(h);
	dims = malloc(numSpaces * sizeof(*dims));
	mboProdSpaceGetDims(h, numSpaces, dims);

	offset = 0;
	embeddingNonZeros(0, numSpaces, dims, blockSize, simpleOp->numFactors,
			  simpleOp->embeddings, rmin, rmax, offset, nnz);

	free(dims);
}

void simpleTOpSparseMatrix(MboProdSpace h, struct SimpleTOp *simpleOp,
			   MboGlobInd rmin, MboGlobInd rmax, int *i, int *j,
			   struct MboAmplitude *a, int *numInserted)
{
	MboGlobInd blockSize;
	MboLocInd *dims;
	int numSpaces;
	struct MboAmplitude alpha;
	
	gatherAllEmbeddings(&simpleOp->numFactors, &simpleOp->embeddings);
	sortEmbeddings(simpleOp->numFactors, simpleOp->embeddings);

	blockSize = mboProdSpaceDim(h);
	numSpaces = mboProdSpaceSize(h);
	dims = malloc(numSpaces * sizeof(*dims));
	mboProdSpaceGetDims(h, numSpaces, dims);

	alpha.re = 1;
	alpha.im = 0;
	embeddingSparseMatrix(0, numSpaces, dims, blockSize, alpha,
			      simpleOp->numFactors, simpleOp->embeddings, rmin,
			      rmax, i, j, a, 0, 0, numInserted);

	free(dims);
}

void simpleTOpDiagonal(MboProdSpace h, struct SimpleTOp *simpleOp,
		       MboGlobInd rmin, MboGlobInd rmax,
		       struct MboAmplitude *diag)
{
	struct MboAmplitude alpha;
	MboGlobInd blockSize;
	int numSpaces;
	MboLocInd *dims;

	gatherAllEmbeddings(&simpleOp->numFactors, &simpleOp->embeddings);
	sortEmbeddings(simpleOp->numFactors, simpleOp->embeddings);

	blockSize = mboProdSpaceDim(h);
	numSpaces = mboProdSpaceSize(h);
	dims = malloc(numSpaces * sizeof(*dims));
	mboProdSpaceGetDims(h, numSpaces, dims);

	alpha.re = 1;
	alpha.im = 0;
	embeddingDiagonal(0, numSpaces, dims, blockSize, alpha,
			  simpleOp->numFactors, simpleOp->embeddings, rmin,
			  rmax, diag, 0, 0);
	free(dims);
}

void simpleTOpDeleteDiagonal(struct SimpleTOp *simpleOp)
{
	MboElemOp null;
	int i;
	if (simpleOp->numFactors == 0) {
		simpleOp->numFactors = 1;
		simpleOp->embeddings = realloc(simpleOp->embeddings,
					       sizeof(*simpleOp->embeddings));
		simpleOp->embeddings[0].i = 0;
		mboElemOpCreate(&null);
		simpleOp->embeddings[0].op = null;
	} else {
		gatherAllEmbeddings(&simpleOp->numFactors,
				    &simpleOp->embeddings);
		sortEmbeddings(simpleOp->numFactors, simpleOp->embeddings);
		for (i = 0; i < simpleOp->numFactors; ++i) {
			embeddingDeleteDiagonal(simpleOp->embeddings + i);
		}
	}
}
