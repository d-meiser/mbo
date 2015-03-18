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
/*
\par
Illustrates the parallel application of an MBO operator.
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "Mbo.h"
#include <config.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

static const int numSpins = 8;
static const int numIters = 1;

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

static MboNumOp buildJx(MboProdSpace h)
{
	MboElemOp sp, sm;
	MboTensorOp Jx;
	MboNumOp Jx_compiled;
	struct MboAmplitude pointFive;
	int i;
	MBO_STATUS err;

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
	err = mboNumOpCompile(Jx, &Jx_compiled);
	assert(err == MBO_SUCCESS);
	mboTensorOpDestroy(&Jx);
	return Jx_compiled;
}

int main()
{
	MboNumOp Jx;
	MboProdSpace hTot;
	struct MboAmplitude one, zero, *x, *y, *yomp;
	int i;
	MboGlobInd n, chunk, chunkSize, numChunks, dim;
	clock_t tstart, tend;
	double deltat, difference;
	int numThreads;
	MboNumSubMatrix *chunks;
	MBO_STATUS err;

	hTot = buildSpace(numSpins);
	dim = mboProdSpaceDim(hTot);
	Jx = buildJx(hTot);

	one.re = 1.0;
	one.im = 0.0;
	zero.re = 0.0;
	zero.im = 0.0;
	x = malloc(dim * sizeof(*x));
	y = malloc(dim * sizeof(*x));
	yomp = malloc(dim * sizeof(*x));
	for (n = 0; n < dim; ++n) {
		x[n].re = (double)rand() / RAND_MAX;
		x[n].im = (double)rand() / RAND_MAX;
	}

	tstart = clock();
	for (i = 0; i < numIters; ++i) {
		mboNumOpMatVec(one, Jx, x, zero, y);
	}
	tend = clock();
	deltat = (double)(tend - tstart) / (double)CLOCKS_PER_SEC;
	printf("Serial: %2.3lf s\n", deltat);

#ifdef HAVE_OPENMP
	numThreads = omp_get_max_threads();
#else
	numThreads = 1;
	printf("here");
#endif
	printf("Using %d threads.\n", numThreads);
	chunkSize = dim / numThreads;
	printf("Chunk Size: %lld\n", chunkSize);
	numChunks = dim / chunkSize;
	printf("Number of Chunks: %lld\n", numChunks);
	if (numChunks * chunkSize < dim) ++numChunks;
	chunks = (MboNumSubMatrix*)malloc(numChunks * sizeof(*chunks));
	for (chunk = 0; chunk < numChunks; ++chunk) {
		err = mboNumSubMatrixCreate(Jx, chunk * chunkSize,
			(chunk + 1) * chunkSize, 0, dim, &chunks[chunk]);
		assert(err == MBO_SUCCESS);
	}
	tstart = clock();
	for (i = 0; i < numIters; ++i) {
#ifdef HAVE_OPENMP
#pragma omp parallel for private(chunk)
#endif
		for (chunk = 0; chunk < numChunks; ++chunk) {
			mboNumSubMatrixMatVec(one, chunks[chunk], x, zero,
					      yomp + chunk * chunkSize);
		}
	}
	tend = clock();
	deltat = (double)(tend - tstart) / (double)CLOCKS_PER_SEC / numThreads;
	printf("Parallel: %2.3lf s\n", deltat);

	difference = 0;
	for (n = 0; n < dim; ++n) {
		difference += (yomp[n].re - y[n].re) * (yomp[n].re - y[n].re) +
			      (yomp[n].im - y[n].im) * (yomp[n].im - y[n].im);
	}
	difference = sqrt(difference);

	printf("Error: %lf\n", difference);

	free(x);
	free(y);
	free(yomp);
	for (chunk = 0; chunk < numChunks; ++chunk) {
		mboNumSubMatrixDestroy(&chunks[chunk]);
	}
	free(chunks);
	mboNumOpDestroy(&Jx);
	mboProdSpaceDestroy(&hTot);

	if (difference < 1.0e-12) {
		return 0;
	} else {
		return 1;
	}

}
