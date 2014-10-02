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

#include <stdlib.h>
#include <MboTensorOpPrivate.h>
#include <SimpleTOp.h>
#include <MboNumOp.h>

struct MboNumOp_t
{
  MboProdSpace space;
  int numTerms;
  struct SimpleTOp *sum;
};

MboNumOp mboNumOpCompile(MboTensorOp op)
{
  struct SimpleTOp *terms;
  int i;
  MboNumOp numOp;

  numOp = malloc(sizeof(*numOp));
  numOp->space = mboProdSpaceCopy(mboTensorOpGetSpace(op)); 
  numOp->numTerms = mboTensorOpGetNumTerms(op);
  numOp->sum = malloc(numOp->numTerms * sizeof(*numOp->sum));
  terms = mboTensorOpGetSimpleTOps(op);
  for (i = 0; i < numOp->numTerms; ++i) {
    numOp->sum[i].numFactors = 0;
    numOp->sum[i].embeddings = 0;
    copySimpleTOp(numOp->sum + i, terms + i);
  }
  return numOp;
}

void mboNumOpDestroy(MboNumOp *op)
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

MboProdSpace mboNumOpGetSpace(MboNumOp op)
{
  return op->space;
}

static void zeroArray(MboGlobInd n, struct MboAmplitude *array)
{
	MboGlobInd i;
	for (i = 0; i < n; ++i) {
		array[i].re = 0;
		array[i].im = 0;
	}
}

void mboNumOpDenseMatrix(MboNumOp a, struct MboAmplitude *mat)
{
	int i;
	MboGlobInd dim;
	dim = mboProdSpaceDim(a->space);

	zeroArray(dim * dim, mat);
	for (i = 0; i < a->numTerms; ++i) {
		simpleTOpDenseMatrix(a->space, a->sum + i, mat);
	}
}

void mboNumOpRowOffsets(MboNumOp op, MboGlobInd rmin, MboGlobInd rmax,
			   int *i)
{
	int j, r;
	for (r = 0; r < rmax - rmin; ++r) {
		i[r] = 0;
	}
	for (j = 0; j < op->numTerms; ++j) {
		simpleTOpGetNonZerosPerRow(op->space, op->sum + j, rmin, rmax,
					   i + 1);
	}
	i[0] = 0;
	for (j = 0; j < rmax - rmin; ++j) {
		i[j + 1] += i[j];
	}
}

void mboNumOpSparseMatrix(MboNumOp op, MboGlobInd rmin, MboGlobInd rmax,
			     int *i, int *j, struct MboAmplitude *a)
{
	int *numInserted, r, s;

	if (rmax <= rmin) return;
	numInserted = malloc((rmax - rmin) * sizeof(*numInserted));
	for (r = 0; r < rmax - rmin; ++r) {
		numInserted[r] = 0;
	}
	for (s = 0; s < op->numTerms; ++s) {
		simpleTOpSparseMatrix(op->space, op->sum + s, rmin, rmax, i, j,
				a, numInserted);
	}
	free(numInserted);
}

MBO_STATUS mboNumOpMatVec(struct MboAmplitude alpha, MboNumOp a,
			     struct MboAmplitude *x, struct MboAmplitude beta,
			     struct MboAmplitude *y, MboGlobInd rmin,
			     MboGlobInd rmax)
{
	int i;
  MboGlobInd r;
	MBO_STATUS err;
	struct MboAmplitude tmp;

  for (r = rmin; r < rmax; ++r) {
    tmp.re = beta.re * y[r].re - beta.im * y[r].im;
    tmp.im = beta.re * y[r].im + beta.im * y[r].re;
    y[r].re = tmp.re;
    y[r].im = tmp.im;
  }

	for (i = 0; i < a->numTerms; ++i) {
		err = applySimpleTOp(a->space, alpha, a->sum + i, x, y, rmin, rmax);
		if (err != MBO_SUCCESS) return err;
	}
	return MBO_SUCCESS;
}

double mboNumOpFlops(MboNumOp a)
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
