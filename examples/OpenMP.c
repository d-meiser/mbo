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
 * \example OpenMP.c
 * Illustrates the parallel application of an MBO operator.
 * */
#include <stdio.h>
#include <stdlib.h>
#include "Mbo.h"

static const int numSpins = 12;
static const int numIters = 10;

int main()
{
	MboElemOp sp, sm;
	MboTensorOp Jx;
	MboProdSpace h1, hTot;
  struct MboAmplitude one, zero, pointFive, *x, *y;
  int i;
  MboGlobInd n;

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

  x = malloc(mboProdSpaceDim(hTot) * sizeof(*x));
  y = malloc(mboProdSpaceDim(hTot) * sizeof(*x));
  for (n = 0; n < mboProdSpaceDim(hTot); ++n) {
    x[n].re = (double)rand() / RAND_MAX;
    x[n].im = (double)rand() / RAND_MAX;
  }
  for (i = 0; i < numIters; ++i) {
    mboTensorOpMatVec(one, Jx, x, zero, y, 0, mboProdSpaceDim(hTot));
  }

	/* Deallocate remaining resources */
  free(x);
  free(y);
	mboTensorOpDestroy(&Jx);
	mboProdSpaceDestroy(&hTot);
}
