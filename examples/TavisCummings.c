#include <Mbo.h>
#include <stdio.h>

static const int nAtoms = 10;
static const int nPhotons = 30;
static const double omegaC = 2.345e0;
static const int numIter = 0;

static double omega(int i);

int main()
{
	MboProdSpace hSingleAtom, hAtoms, hField, hTot;
	MboElemOp sm, sp, sz, a, ad, numberOp;
	MboTensorOp inhomogeneousJz, jPlus, jMinus, idField, idAtoms, A, Ad, N,
	    H;
	struct MboAmplitude tmp;
	MboVec x, y;
	int i;

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
	for (int i = 0; i < nAtoms; ++i) {
		tmp.re = omega(i);
		tmp.im = 0;
		mboTensorOpAddScaledTo(&tmp, sz, i, inhomogeneousJz);
	}
	mboTensorOpNull(hAtoms, &jMinus);
	for (int i = 0; i < nAtoms; ++i) {
		mboTensorOpAddTo(sm, i, jMinus);
	}
	mboTensorOpNull(hAtoms, &jPlus);
	for (int i = 0; i < nAtoms; ++i) {
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
	mboTensorOpNull(hTot, &H);
	mboTensorOpKron(idField, inhomogeneousJz, &H);
	mboTensorOpKron(N, idAtoms, &H);
	mboTensorOpKron(A, jPlus, &H);
	mboTensorOpKron(Ad, jMinus, &H);

	printf("Dimension of total Hilbert Space: %lld\n",
	       mboProdSpaceDim(hTot));

	mboVecCreate(mboProdSpaceDim(hTot), &x);
	mboVecCreate(mboProdSpaceDim(hTot), &y);

	tmp.re = 1.0;
	tmp.im = 0.0;
	mboTensorOpMatVec(&tmp, H, x, &tmp, y);
	for (i = 0; i < numIter; ++i) {
		mboTensorOpMatVec(&tmp, H, x, &tmp, y);
	}

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

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}

