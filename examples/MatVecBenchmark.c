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
These benchmarks compare the performance of matrix vector multiplies
using #MboNumOp to assembled sparse matrix vector multiplies.  The
assembled matrices correspond to regular compressed sparse row (CSR)
matrices.  This matrix format is sometimes referred to as Yale formate.


<h3>Performance model for assembled operator application</h3>

\par
Operator application for assembled matrices is typically memory
bandwidth limited.  For each complex multiply-add we need to load at
least \f$(8+16)\f$ Byte.  The first \f$8\f$ Byte are for the column
index of the non-zero.  The second \f$16\f$ Byte are the complex number
for the non-zero entry in the matrix.  Note that this estimate assumes
that the right hand side values are caches which is not always the case.
24B per multiply add is thus an optimistic estimate.  It corresponds to
an arithmetic intensity of \f$1/3\f$ floating point operations per Byte.
For a stream bandwidth of \f$ {\rm BW} \f$ we get \f$ 1/3 ({\rm BW} /
({\rm GB}/{\rm s})) {\rm GFLOP}/{\rm s} \f$.


<h3>Performance model for #MboNumOp</h3>

\par
The implementation of the matrix vector multiply inside of #MboNumOp is
going to evolve.  Currently, this method is essentially streaming both
source and destination vector for each simple operator.  A simple
operator in this context is a tensor product composed entirely of
identities and individual sparse matrix entries.  This means that the
#MboNumOp matrix vector multiply is largely memory bandwidth limited.
In its current implementation we have an arithmetic intensity of
\f$1/4\f$.  The fact that this implementation is still slightly faster
than the assembled matrix vector multiply despite being even more memory
bandwidth bound is probably due to the more regular data access pattern.
We avoid the indirect data accesses in the assembled matrix vector
multiply.


\par
Before we start we need to bring the Mbo library into scop by including
several headers.  We need the main Mbo.h header as well as MboVec.h
*/
#include <Mbo.h>
#include <MboVec.h>
/*
\par
In addition we need the following system headers
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

/*
\par
The following two functions are used to construct the #MboNumOp's needed
in this benchmark.  They are defined further down.
*/
static MboNumOp buildTavisCummingsHamiltonian(int nAtoms, int nPhotons);
static MboNumOp buildJx(int nAtoms);

/*
\par
The following function builds all operators we're using for the
comparison.  Currently we're using the many particle angular momentum
operator \f$J_x\f$ as well as the Tavis-Cummings Hamiltonian.  `numOps`
specifies the total number of operators that have been allocated in the
`ops` array.
*/
static void buildAllOps(int numOps, MboNumOp* ops)
{
	int i;

	static const int numPtclsJxMin = 8;
	static const int numPtclsJxMax = 18;
	for (i = numPtclsJxMin; i < numPtclsJxMax; ++i) {
		if (numOps <= 0) return;
		*ops = buildJx(i);
		++ops;
		--numOps;
	}
	for (i = 0; i < 10; ++i) {
		if (numOps <= 0) return;
		*ops = buildTavisCummingsHamiltonian(10 + i, 10);
		++ops;
		--numOps;
	}
}

static void runBenchmark(MboNumOp op);

/*
\par
The benchmarks are run by simply creating the operators and then
invoking runBenchmark on each one of them.
*/
int main() {
	MboNumOp *ops;
	int numOps = 20, i;

	ops = (MboNumOp*)malloc(numOps * sizeof(*ops));
	buildAllOps(numOps, ops);

	for (i = 0; i < numOps; ++i) {
		runBenchmark(ops[i]);
	}

/*
Here we release all resources associated with the operators.
*/
	for (i = 0; i < numOps; ++i) {
		mboNumOpDestroy(&ops[i]);
	}
	free(ops);
}

/*
\par
The CSR data structure defined in the following code is used to
represent the assembled matrices.  `dim` is the dimension of the matrix;
`i` is an array with the offsets for the different rows; `j` are the
column indices of the non-zero entries; and `a` are the non-zero values.
*/
struct CSR {
	MboGlobInd dim;
	MboGlobInd *i;
	MboGlobInd *j;
	struct MboAmplitude *a;
};

/*
\par
The following function initializes a CSR structure in a well defined
state.
*/
static struct CSR createCSR()
{
	struct CSR csr;
	csr.i = 0;
	csr.j = 0;
	csr.a = 0;
	return csr;
}

/*
\par
The following function deallocates all memory in a CSR matrix.
*/
static void destroyCSR(struct CSR *csr)
{
	free(csr->i);
	free(csr->j);
	free(csr->a);
}

/*
\par
`getCSRMatrix` creates the CSR matrix corresponding to a #MboNumOp.
This function uses #mboNumOpRowOffsets and #mboNumOpSparseMatrix to
obtain the sparse matrix data.  `csr` should have been previously
created with `createCSR`.
*/
static void getCSRMatrix(MboNumOp op, struct CSR *csr)
{
	MboProdSpace h;
	MboGlobInd dim;

	h = mboNumOpGetSpace(op);
	dim = mboProdSpaceDim(h);
	csr->dim = dim;
	csr->i = malloc((dim + 1) * sizeof(*csr->i));
	mboNumOpRowOffsets(op, 0, dim, csr->i);
	csr->j = malloc(csr->i[dim] * sizeof(*csr->j));
	csr->a = malloc(csr->i[dim] * sizeof(*csr->a));
	mboNumOpSparseMatrix(op, 0, dim, csr->i, csr->j, csr->a);
}

/*
\par
`BenchmarkData` is a simple container for the results from a benchmark run.  It
is used to record the execution times of several runs, storage requirements, and
the number of floating point operations.
*/
struct BenchmarkData {
	int numIters;
	double *times;
	double storage;
	double numFlops;
};

/*
\par
A `struct BenchmarkData` can be initialized with the following function
*/
static struct BenchmarkData createBenchmarkData(int numIters)
{
	struct BenchmarkData bm;
	bm.numIters = numIters;
	bm.times = malloc(numIters * sizeof(*bm.times));
	bm.storage = 0;
	bm.numFlops = 0;
	return bm;
}

/*
and it is destroyed with
*/
static void destroyBenchmarkData(struct BenchmarkData *bm)
{
	free(bm->times);
}

/*
\par
A CSR matrix can be applied to a vector by means of the following function.
*/
static void applyCSRMatrix(struct MboAmplitude alpha, struct CSR *csr,
			   struct MboAmplitude *x, struct MboAmplitude beta,
			   struct MboAmplitude *y)
{
	MboGlobInd row, col, i;
	struct MboAmplitude tmp, tmpp;

	for (row = 0; row < csr->dim; ++row) {
		tmp.re = beta.re * y[row].re - beta.im * y[row].im;
		tmp.im = beta.re * y[row].im + beta.im * y[row].re;
		for (i = csr->i[row]; i < csr->i[row + 1]; ++i) {
			col = csr->j[i];
			tmpp.re = csr->a[i].re * alpha.re - csr->a[i].im *
				alpha.im;
			tmpp.im = csr->a[i].re * alpha.im + csr->a[i].im *
				alpha.re;
			tmp.re += x[col].re * tmpp.re - x[col].im * tmpp.im;
			tmp.im += x[col].re * tmpp.im + x[col].im * tmpp.re;
		}
		y[row].re = tmp.re;
		y[row].im = tmp.im;
	}
}

/*
\par
This is a simple driver for timing the execution time of the matrix vector
multiplication using assembled matrices.  This driver simply invokes
`applyCSRMatrix` repeatedly and records the execution time.
*/
static struct BenchmarkData runAssembledBenchmark(MboNumOp op, int numIters)
{
	struct CSR csr = createCSR();
	struct BenchmarkData bm = createBenchmarkData(numIters);
	MboVec x, y;
	struct MboAmplitude *xarr, *yarr, alpha, beta;
	int i;
	clock_t tstart, tend;

	getCSRMatrix(op, &csr);
	bm.numFlops = csr.i[csr.dim] * 8.0;

	mboVecCreate(csr.dim, &x);
	mboVecCreate(csr.dim, &y);
	mboVecGetViewR(x, &xarr);
	mboVecGetViewRW(y, &yarr);
	alpha.re = 2.7;
	alpha.im = 0.8;
	beta.re = 2.1;
	beta.im = 0.9;
	for (i = 0; i < numIters; ++i) {
		tstart = clock();
		applyCSRMatrix(alpha, &csr, xarr, beta, yarr);
		tend = clock();
		bm.times[i] = (double)(tend - tstart) / CLOCKS_PER_SEC;
	}

	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);
	
	mboVecDestroy(&x);
	mboVecDestroy(&y);
	destroyCSR(&csr);
	return bm;
}

/*
The driver for the #MboNumOp's is completely analogous:
*/
static struct BenchmarkData runMboBenchmark(MboNumOp op, int numIters)
{
	clock_t tstart, tend;
	MboVec x, y;
	struct MboAmplitude *xarr, *yarr, alpha, beta;
	struct BenchmarkData bm = createBenchmarkData(numIters);
	MboGlobInd dim;
	int i;

	dim = mboProdSpaceDim(mboNumOpGetSpace(op));
	bm.numFlops = mboNumOpFlops(op);

	mboVecCreate(dim, &x);
	mboVecCreate(dim, &y);
	mboVecGetViewR(x, &xarr);
	mboVecGetViewRW(y, &yarr);
	alpha.re = 2.7;
	alpha.im = 0.8;
	beta.re = 2.1;
	beta.im = 0.9;
	for (i = 0; i < numIters; ++i) {
		tstart = clock();
		mboNumOpMatVec(alpha, op, xarr, beta, yarr);
		tend = clock();
		bm.times[i] = (double)(tend - tstart) / CLOCKS_PER_SEC;
	}

	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);
	
	mboVecDestroy(&x);
	mboVecDestroy(&y);
	return bm;
}

/*
\par
The average of several numbers is computed as follows.  This can be used to get
average run times for the different benchmarks.
*/
static double average(const double *x, int n)
{
	double avg = 0.0;
	int i;

	for (i = 0; i < n; ++i) {
		avg += x[i];
	}

	return avg / (double)n;
}

/*
\par
A pair of assembled/unassembled benchmarks is printed to the screen as follows
*/
static void printBenchmarkReport(struct BenchmarkData *assembledBM,
				 struct BenchmarkData *mboBenchmark, int valid)
{
	double meanTimeAssembled, meanTimeMbo;
	double gflopsAssembled, gflopsMbo;

	meanTimeAssembled = average(assembledBM->times, assembledBM->numIters);
	meanTimeMbo = average(mboBenchmark->times, mboBenchmark->numIters);


	gflopsAssembled = assembledBM->numFlops / meanTimeAssembled / 1.0e9;
	gflopsMbo = assembledBM->numFlops / meanTimeMbo / 1.0e9;
	printf("%1.2lf (%1.2Es)   %1.2lf (%1.2Es)   %1.2lfx ", gflopsMbo,
	       meanTimeMbo, gflopsAssembled, meanTimeMbo,
	       meanTimeAssembled / meanTimeMbo);
	if (valid) {
		printf("[PASSED]\n");
	} else {
		printf("[FAILED]\n");
	}
}

/*
\par
To compute the error we use the L2 norm of the difference between two vectors:
*/
static double distance(MboVec x, MboVec y, MboGlobInd d)
{
	MboGlobInd i;
	double dist = 0;
	struct MboAmplitude diff, *xarr, *yarr;

	mboVecGetViewR(x, &xarr);
	mboVecGetViewR(y, &yarr);
	for (i = 0; i < d; ++i) {
		diff.re = xarr[i].re - yarr[i].re;
		diff.im = xarr[i].im - yarr[i].im;
		dist += diff.re * diff.re + diff.im * diff.im;
	}
	dist = sqrt(dist);
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);
	return dist;
}

/*
\par
The `verify` method compares the results obtained from the assembled operators
and the unassembled operators in the L2 norm.  If the error is "small enough" we
return a nonzero value to indicate success.
*/
static int verify(MboNumOp op)
{
	int valid = 0;
	MboVec x, yMboNumOp, yAssembled;
	struct MboAmplitude *xarr, *yarr, alpha, beta;
	MboGlobInd dim;
	struct CSR csr = createCSR();
	double error;

	dim = mboProdSpaceDim(mboNumOpGetSpace(op));
	alpha.re = 2.7;
	alpha.im = 0.8;
	beta.re = 0;
	beta.im = 0;
	mboVecCreate(dim, &x);
	mboVecSetRandom(x);

/*
Compute result with MboNumOp
*/
	mboVecCreate(dim, &yMboNumOp);
	mboVecSet(&beta, yMboNumOp);
	mboVecGetViewR(x, &xarr);
	mboVecGetViewRW(yMboNumOp, &yarr);
	mboNumOpMatVec(alpha, op, xarr, beta, yarr);
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(yMboNumOp, &yarr);

/*
Compute result with assembled operator
*/
	mboVecCreate(dim, &yAssembled);
	mboVecSet(&beta, yAssembled);
	mboVecGetViewR(x, &xarr);
	mboVecGetViewRW(yAssembled, &yarr);
	getCSRMatrix(op, &csr);
	applyCSRMatrix(alpha, &csr, xarr, beta, yarr);
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(yAssembled, &yarr);

	error = distance(yMboNumOp, yAssembled, dim) / sqrt((double)dim);

	mboVecDestroy(&x);
	mboVecDestroy(&yMboNumOp);
	mboVecDestroy(&yAssembled);
	destroyCSR(&csr);

	static const double TOLERANCE = 1.0e-10;
	if (error > TOLERANCE) {
		valid = 0;
	} else {
		valid = 1;
	}

	return valid;
}

/*
\par
Running of the benchmarks, verification, and reporting for a given operator is
orchestrated as follows:
*/
void runBenchmark(MboNumOp op)
{
	static const int numIters = 2;
	struct BenchmarkData assembledBenchmark =
	    runAssembledBenchmark(op, numIters);
	struct BenchmarkData mboBenchmark = runMboBenchmark(op, numIters);
	int operatorValid = verify(op);
	printBenchmarkReport(&assembledBenchmark, &mboBenchmark, operatorValid);
	destroyBenchmarkData(&assembledBenchmark);
	destroyBenchmarkData(&mboBenchmark);
}

/*
\par
Resonance frequency of ith atom used in the construction of the Tavis-Cummings
Hamiltonian below.  This is an entirely arbitrary choice
*/
static double omega(int i)
{
	return 1.3 * i;
}

/*
\par
Finally, the last two functions are used to construct the operators in the
benchmark.
*/
MboNumOp buildTavisCummingsHamiltonian(int nAtoms, int nPhotons)
{
	MboProdSpace hSingleAtom, hAtoms, hField, hTot;
	MboElemOp sm, sp, sz, a, ad, numberOp;
	MboTensorOp inhomogeneousJz, jPlus, jMinus, idField, idAtoms, A, Ad, N,
	    H, *factors;
	MboNumOp H_compiled;
	struct MboAmplitude tmp;
	int i;
	MBO_STATUS err;

/*
build various Hilbert spaces
*/
	hSingleAtom = mboProdSpaceCreate(2);
	hAtoms = mboProdSpaceCreate(0);
	for (i = 0; i < nAtoms; ++i) {
		mboProdSpaceMul(hSingleAtom, &hAtoms);
	}
	hField = mboProdSpaceCreate(nPhotons + 1);
	hTot = mboProdSpaceCreate(0);
	mboProdSpaceMul(hAtoms, &hTot);
	mboProdSpaceMul(hField, &hTot);

/*
build atomic operators
*/
	sm = mboSigmaMinus();
	sp = mboSigmaPlus();
	sz = mboSigmaZ();
	mboTensorOpNull(hAtoms, &inhomogeneousJz);
	for (i = 0; i < nAtoms; ++i) {
		tmp.re = omega(i);
		tmp.im = 0;
		mboTensorOpAddScaledTo(&tmp, sz, i, inhomogeneousJz);
	}
	mboTensorOpNull(hAtoms, &jMinus);
	for (i = 0; i < nAtoms; ++i) {
		mboTensorOpAddTo(sm, i, jMinus);
	}
	mboTensorOpNull(hAtoms, &jPlus);
	for (i = 0; i < nAtoms; ++i) {
		mboTensorOpAddTo(sp, i, jPlus);
	}
	mboTensorOpIdentity(hAtoms, &idAtoms);

/*
build field operators
*/
	a = mboAnnihilationOp(nPhotons + 1);
	ad = mboCreationOp(nPhotons + 1);
	numberOp = mboNumOp(nPhotons + 1);
	mboTensorOpNull(hField, &A);
	mboTensorOpAddTo(a, 0, A);
	mboTensorOpNull(hField, &Ad);
	mboTensorOpAddTo(ad, 0, Ad);
	mboTensorOpNull(hField, &N);
/*
The frequency of the cavity
*/
	tmp.re = 2.7;
	tmp.im = 0;
	mboTensorOpAddScaledTo(&tmp, numberOp, 0, N);
	mboTensorOpIdentity(hField, &idField);

/*
Assemble Hamiltonian
*/
	factors = malloc(2 * sizeof(*factors));
	mboTensorOpNull(hTot, &H);
	factors[0] = idField;
	factors[1] = inhomogeneousJz;
	mboTensorOpKron(2, factors, &H);
	factors[0] = N;
	factors[1] = idAtoms;
	mboTensorOpKron(2, factors, &H);
	factors[0] = A;
	factors[1] = jPlus;
	mboTensorOpKron(2, factors, &H);
	factors[0] = Ad;
	factors[1] = jMinus;
	mboTensorOpKron(2, factors, &H);
	free(factors);

	err = mboNumOpCompile(H, &H_compiled);
	assert(err == MBO_SUCCESS);

/*
Release all resources
*/
	mboElemOpDestroy(&sm);
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sz);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&ad);
	mboElemOpDestroy(&numberOp);
	mboProdSpaceDestroy(&hSingleAtom);
	mboProdSpaceDestroy(&hAtoms);
	mboProdSpaceDestroy(&hField);
	mboProdSpaceDestroy(&hTot);
	mboTensorOpDestroy(&inhomogeneousJz);
	mboTensorOpDestroy(&jMinus);
	mboTensorOpDestroy(&jPlus);
	mboTensorOpDestroy(&idAtoms);
	mboTensorOpDestroy(&A);
	mboTensorOpDestroy(&Ad);
	mboTensorOpDestroy(&N);
	mboTensorOpDestroy(&idField);
	mboTensorOpDestroy(&H);

	return H_compiled;
}

MboNumOp buildJx(int nAtoms)
{
	MboProdSpace hSingleAtom, hAtoms;
	MboElemOp sm, sp;
	MboTensorOp Jx;
	MboNumOp Jx_compiled;
	struct MboAmplitude tmp;
	int i;
	MBO_STATUS err;

/*
build various Hilbert spaces
*/
	hSingleAtom = mboProdSpaceCreate(2);
	hAtoms = mboProdSpaceCreate(0);
	for (i = 0; i < nAtoms; ++i) {
		mboProdSpaceMul(hSingleAtom, &hAtoms);
	}

/* 
build atomic operators
*/
	sm = mboSigmaMinus();
	sp = mboSigmaPlus();
	mboTensorOpNull(hAtoms, &Jx);
	for (i = 0; i < nAtoms; ++i) {
		tmp.re = 0.5;
		tmp.im = 0;
		mboTensorOpAddScaledTo(&tmp, sm, i, Jx);
		mboTensorOpAddScaledTo(&tmp, sp, i, Jx);
	}

	err = mboNumOpCompile(Jx, &Jx_compiled);
	assert(err == MBO_SUCCESS);

/*
Release all resources
*/
	mboElemOpDestroy(&sm);
	mboElemOpDestroy(&sp);
	mboProdSpaceDestroy(&hSingleAtom);
	mboProdSpaceDestroy(&hAtoms);
	mboTensorOpDestroy(&Jx);

	return Jx_compiled;
}
