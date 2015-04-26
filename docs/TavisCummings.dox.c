/**
@defgroup TavisCummings TavisCummings
@ingroup mbo_examples
@{
\par
Illustrates the construction of the Tavis-Cummings Hamiltonian using
MBO.

\code
#include <Mbo.h>
#include <MboVec.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static const int nAtoms = 10;
static const int nPhotons = 30;
static const double omegaC = 2.345e0;
static const int numIter = 2;

static double omega(int i);

int main()
{
	MboProdSpace hSingleAtom, hAtoms, hField, hTot;
	MboElemOp sm, sp, sz, a, ad, numberOp;
	MboTensorOp inhomogeneousJz, jPlus, jMinus, idField, idAtoms, A, Ad, N,
	    H, *factors;
	MboNumOp H_compiled;
	struct MboAmplitude tmp, *xarr, *yarr;
	MboVec x, y;
	int i;
	MBO_STATUS err;

	/* build various Hilbert spaces */
	hSingleAtom = mboProdSpaceCreate(2);
	hAtoms = mboProdSpaceCreate(0);
	for (i = 0; i < nAtoms; ++i) {
		mboProdSpaceMul(hSingleAtom, &hAtoms);
	}
	hField = mboProdSpaceCreate(nPhotons + 1);
	hTot = mboProdSpaceCreate(0);
	mboProdSpaceMul(hAtoms, &hTot);
	mboProdSpaceMul(hField, &hTot);

	/* build atomic operators */
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

	/* build field operators */
	a = mboAnnihilationOp(nPhotons + 1);
	ad = mboCreationOp(nPhotons + 1);
	numberOp = mboNumOp(nPhotons + 1);
	mboTensorOpNull(hField, &A);
	mboTensorOpAddTo(a, 0, A);
	mboTensorOpNull(hField, &Ad);
	mboTensorOpAddTo(ad, 0, Ad);
	mboTensorOpNull(hField, &N);
	/* The frequency of the cavity */
	tmp.re = omegaC;
	tmp.im = 0;
	mboTensorOpAddScaledTo(&tmp, numberOp, 0, N);
	mboTensorOpIdentity(hField, &idField);

	/* Assemble Hamiltonian */
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

	printf("Dimension of total Hilbert Space: %lld\n",
	       mboProdSpaceDim(hTot));
	printf("GFLOPs per operator application (estimated):  %lf\n",
	       mboNumOpFlops(H_compiled) / (double)(1 << 30));
	printf("Total GFLOPs (estimated):  %lf\n",
	       (1 + numIter) * mboNumOpFlops(H_compiled) / (double)(1 << 30));

	mboVecCreate(mboProdSpaceDim(hTot), &x);
	mboVecCreate(mboProdSpaceDim(hTot), &y);

	tmp.re = 1.0;
	tmp.im = 0.0;
	mboVecGetViewR(x, &xarr);
	mboVecGetViewR(y, &yarr);
	mboNumOpMatVec(tmp, H_compiled, xarr, tmp, yarr);
	for (i = 0; i < numIter; ++i) {
		mboNumOpMatVec(tmp, H_compiled, xarr, tmp, yarr);
	}
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);

	/* Release all resources */
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
	mboVecDestroy(&x);
	mboVecDestroy(&y);
  mboNumOpDestroy(&H_compiled);

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}

\endcode

<h2>Plain source code</h2>


\code

#include <Mbo.h>
#include <MboVec.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static const int nAtoms = 10;
static const int nPhotons = 30;
static const double omegaC = 2.345e0;
static const int numIter = 2;

static double omega(int i);

int main()
{
	MboProdSpace hSingleAtom, hAtoms, hField, hTot;
	MboElemOp sm, sp, sz, a, ad, numberOp;
	MboTensorOp inhomogeneousJz, jPlus, jMinus, idField, idAtoms, A, Ad, N,
	    H, *factors;
	MboNumOp H_compiled;
	struct MboAmplitude tmp, *xarr, *yarr;
	MboVec x, y;
	int i;
	MBO_STATUS err;

	/* build various Hilbert spaces */
	hSingleAtom = mboProdSpaceCreate(2);
	hAtoms = mboProdSpaceCreate(0);
	for (i = 0; i < nAtoms; ++i) {
		mboProdSpaceMul(hSingleAtom, &hAtoms);
	}
	hField = mboProdSpaceCreate(nPhotons + 1);
	hTot = mboProdSpaceCreate(0);
	mboProdSpaceMul(hAtoms, &hTot);
	mboProdSpaceMul(hField, &hTot);

	/* build atomic operators */
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

	/* build field operators */
	a = mboAnnihilationOp(nPhotons + 1);
	ad = mboCreationOp(nPhotons + 1);
	numberOp = mboNumOp(nPhotons + 1);
	mboTensorOpNull(hField, &A);
	mboTensorOpAddTo(a, 0, A);
	mboTensorOpNull(hField, &Ad);
	mboTensorOpAddTo(ad, 0, Ad);
	mboTensorOpNull(hField, &N);
	/* The frequency of the cavity */
	tmp.re = omegaC;
	tmp.im = 0;
	mboTensorOpAddScaledTo(&tmp, numberOp, 0, N);
	mboTensorOpIdentity(hField, &idField);

	/* Assemble Hamiltonian */
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

	printf("Dimension of total Hilbert Space: %lld\n",
	       mboProdSpaceDim(hTot));
	printf("GFLOPs per operator application (estimated):  %lf\n",
	       mboNumOpFlops(H_compiled) / (double)(1 << 30));
	printf("Total GFLOPs (estimated):  %lf\n",
	       (1 + numIter) * mboNumOpFlops(H_compiled) / (double)(1 << 30));

	mboVecCreate(mboProdSpaceDim(hTot), &x);
	mboVecCreate(mboProdSpaceDim(hTot), &y);

	tmp.re = 1.0;
	tmp.im = 0.0;
	mboVecGetViewR(x, &xarr);
	mboVecGetViewR(y, &yarr);
	mboNumOpMatVec(tmp, H_compiled, xarr, tmp, yarr);
	for (i = 0; i < numIter; ++i) {
		mboNumOpMatVec(tmp, H_compiled, xarr, tmp, yarr);
	}
	mboVecReleaseView(x, &xarr);
	mboVecReleaseView(y, &yarr);

	/* Release all resources */
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
	mboVecDestroy(&x);
	mboVecDestroy(&y);
  mboNumOpDestroy(&H_compiled);

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}


\endcode
@}
*/
