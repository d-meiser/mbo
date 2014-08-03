#include <Mbo.h>
#include <stdio.h>

double omega(int i);

int main()
{
	MboProdSpace hSingleAtom;
	MboProdSpace hAtoms;
	MboProdSpace hField;
	MboProdSpace hTot;
	MboElemOp sm;
	MboElemOp sp;
	MboElemOp sz;
	const int nAtoms = 20;
	const int nPhotons = 30;
	MboTensorOp inhomogeneousJz;
	MboTensorOp jMinus;
  struct MboAmplitude tmp;

	mboElemOpCreate(&sm);
  tmp.re = 1.0;
  tmp.im = 0.0;
	mboElemOpAddTo(0, 1, &tmp, &sm);
	mboElemOpCreate(&sp);
	mboElemOpAddTo(1, 0, &tmp, &sp);
	mboElemOpCreate(&sz);
	mboElemOpAddTo(1, 1, &tmp, &sz);
  tmp.re = -1.0;
	mboElemOpAddTo(0, 0, &tmp, &sz);

	hSingleAtom = mboProdSpaceCreate(2);
	hAtoms = mboProdSpaceCreate(0);
	for (int i = 0; i < nAtoms; ++i) {
		mboProdSpaceMul(hSingleAtom, &hAtoms);
	}
	hField = mboProdSpaceCreate(nPhotons + 1);
	hTot = mboProdSpaceCreate(0);
	mboProdSpaceMul(hAtoms, &hTot);
	mboProdSpaceMul(hField, &hTot);

	mboTensorOpNull(hAtoms, &inhomogeneousJz);
	for (int i = 0; i < nAtoms; ++i) {
    tmp.re = omega(i);
		mboTensorOpAddScaledTo(&tmp, sz, i, inhomogeneousJz);
	}

	mboTensorOpNull(hTot, &jMinus);

	mboElemOpDestroy(&sm);
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sz);
	mboProdSpaceDestroy(&hSingleAtom);
	mboProdSpaceDestroy(&hAtoms);
	mboProdSpaceDestroy(&hField);
	mboProdSpaceDestroy(&hTot);
	mboTensorOpDestroy(&inhomogeneousJz);
	mboTensorOpDestroy(&jMinus);

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}

