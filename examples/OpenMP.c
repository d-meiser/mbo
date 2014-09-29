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
#include <time.h>
#include <math.h>
#include "Mbo.h"

static const int numSpins = 13;
static const int numIters = 10;

static MboProdSpace buildSpace(int n)
{
	MboProdSpace h1, hTot;
	int i;

	h1 = mboProdSpaceCreate(2);
	hTot = mboProdSpaceCreate(0);
	for (i = 0; i < numSpins; ++i) {
		mboProdSpaceMul(h1, &hTot);
	}
	mboProdSpaceDestroy(&h1);
	return hTot;
}

static MboTensorOp buildJx(MboProdSpace h)
{
	MboElemOp sp, sm;
	MboTensorOp Jx;
	struct MboAmplitude pointFive;
	int i;

	sp = mboSigmaPlus();
	sm = mboSigmaMinus();
	pointFive.re = 0.5;
	pointFive.im = 0.0;
	mboTensorOpNull(h, &Jx);
	for (i = 0; i < numSpins; ++i) {
		mboTensorOpAddScaledTo(&pointFive, sm, i, Jx);
		mboTensorOpAddScaledTo(&pointFive, sp, i, Jx);
	}
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sm);
	return Jx;
}

int main()
{
	MboTensorOp Jx;
	MboProdSpace hTot;
	struct MboAmplitude one, zero, *x, *y, *yomp;
	int i;
	MboGlobInd n, chunk, chunkSize, numChunks;
	clock_t tstart, tend;
	double deltat, difference;

	hTot = buildSpace(numSpins);
	Jx = buildJx(hTot);

	one.re = 1.0;
	one.im = 0.0;
	zero.re = 0.0;
	zero.im = 0.0;
	x = malloc(mboProdSpaceDim(hTot) * sizeof(*x));
	y = malloc(mboProdSpaceDim(hTot) * sizeof(*x));
	yomp = malloc(mboProdSpaceDim(hTot) * sizeof(*x));
	for (n = 0; n < mboProdSpaceDim(hTot); ++n) {
		x[n].re = (double)rand() / RAND_MAX;
		x[n].im = (double)rand() / RAND_MAX;
	}


	tstart = clock();
	for (i = 0; i < numIters; ++i) {
		mboTensorOpMatVec(one, Jx, x, zero, y, 0,
				  mboProdSpaceDim(hTot));
	}
	tend = clock();
	deltat = (double)(tend - tstart) / (double)CLOCKS_PER_SEC;
	printf("Serial: %2.3lf s\n", deltat);

	tstart = clock();
	chunkSize = 1024;
	printf("Chunk Size: %lld\n", chunkSize);
	numChunks = mboProdSpaceDim(hTot) / chunkSize;
	printf("Number of Chunks: %lld\n", numChunks);
	if (numChunks * chunkSize < mboProdSpaceDim(hTot)) ++numChunks;
	for (i = 0; i < numIters; ++i) {
		for (chunk = 0; chunk < numChunks; ++chunk) {
			mboTensorOpMatVec(one, Jx, x + chunk * chunkSize, zero,
					  yomp + chunk * chunkSize,
					  chunk * chunkSize,
					  (chunk + 1) * chunkSize);
		}
	}
	tend = clock();
	deltat = (double)(tend - tstart) / (double)CLOCKS_PER_SEC;
	printf("Parallel: %2.3lf s\n", deltat);

	difference = 0;
	for (n = 0; n < mboProdSpaceDim(hTot); ++n) {
		difference += (yomp[n].re - y[n].re) * (yomp[n].re - y[n].re) +
			      (yomp[n].im - y[n].im) * (yomp[n].im - y[n].im);
	}
	difference = sqrt(difference);
	printf("Error: %lf\n", difference);
	if (difference < 1.0e-12) {
		return 0;
	} else {
		return 1;
	}

	free(x);
	free(y);
	free(yomp);
	mboTensorOpDestroy(&Jx);
	mboProdSpaceDestroy(&hTot);
}
