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
/**
\par
Illustrates the conversion of a MBO operator to a dense matrix.
*/
#include <stdio.h>
#include <assert.h>
#include "MboVec.h"
#include "Mbo.h"

static const int numSpins = 4;

int main()
{
	int i, j;
	MboVec x, y, result;
	MboElemOp sp, sm;
	MboTensorOp Jx;
	MboNumOp Jx_compiled;
	MboProdSpace h1, hTot;
	struct MboAmplitude pointFive, one, zero, dotProduct, *yarr, *resultarr;
	MBO_STATUS err;

	/* Build Hilbert space */
	h1 = mboProdSpaceCreate(2);
	hTot = mboProdSpaceCreate(0);
	for (i = 0; i < numSpins; ++i) {
		mboProdSpaceMul(h1, &hTot);
	}
	mboProdSpaceDestroy(&h1);
	
	/* Build Jx */
	sp = mboSigmaPlus();
	sm = mboSigmaMinus();
	pointFive.re = 0.5;
	pointFive.im = 0.0;
	mboTensorOpNull(hTot, &Jx);
	for (i = 0; i < numSpins; ++i) {
		mboTensorOpAddScaledTo(&pointFive, sm, i, Jx);
		mboTensorOpAddScaledTo(&pointFive, sp, i, Jx);
	}
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sm);

	/* Compute matrix */
	one.re = 1.0;
	one.im = 0.0;
	zero.re = 0.0;
	zero.im = 0.0;
	mboVecCreate(mboProdSpaceDim(hTot), &x);
	mboVecCreate(mboProdSpaceDim(hTot), &y);
	mboVecCreate(mboProdSpaceDim(hTot), &result);

	err = mboNumOpCompile(Jx, &Jx_compiled);
	assert(err == MBO_SUCCESS);

	printf("\nMatrix entries:\n");
	for (i = 0; i < mboProdSpaceDim(hTot); ++i) {
		mboVecUnitVector(i, x);
		for (j = 0; j < mboProdSpaceDim(hTot); ++j) {
			mboVecUnitVector(j, y);
			mboVecGetViewR(y, &yarr);
			mboVecGetViewRW(result, &resultarr);
			mboNumOpMatVec(one, Jx_compiled, yarr, zero, resultarr); 
			mboVecReleaseView(y, &yarr);
			mboVecReleaseView(result, &resultarr);
			mboVecDot(x, result, &dotProduct);
			printf("%1.1lf %1.1lf  ", dotProduct.re, dotProduct.im);
		}
		printf("\n");
	}
	printf("\nMatrix pattern:\n");
	for (i = 0; i < mboProdSpaceDim(hTot); ++i) {
		mboVecUnitVector(i, x);
		for (j = 0; j < mboProdSpaceDim(hTot); ++j) {
			mboVecUnitVector(j, y);
			mboVecGetViewR(y, &yarr);
			mboVecGetViewRW(result, &resultarr);
			mboNumOpMatVec(one, Jx_compiled, yarr, zero, resultarr); 
			mboVecReleaseView(y, &yarr);
			mboVecReleaseView(result, &resultarr);
			mboVecDot(x, result, &dotProduct);
			if (dotProduct.re * dotProduct.re +
				dotProduct.im * dotProduct.im >
			    1.0e-6) {
				printf("x");
			} else {
				printf(" ");
			}
		}
		printf("\n");
	}
	mboVecDestroy(&x);
	mboVecDestroy(&y);
	mboVecDestroy(&result);

	/* Deallocate remaining resources */
  mboNumOpDestroy(&Jx_compiled);
	mboTensorOpDestroy(&Jx);
	mboProdSpaceDestroy(&hTot);
}
