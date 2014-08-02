#include <Mbo.h>
#include <stdio.h>

double omega(int i);

int main()
{
	ProdSpace hSingleAtom;
	ProdSpace hAtoms;
	ProdSpace hField;
	ProdSpace hTot;
	ElemOp sm;
	ElemOp sp;
	ElemOp sz;
	const int nAtoms = 20;
	const int nPhotons = 30;
	TensorOp inhomogeneousJz;
	TensorOp jMinus;

	elemOpCreate(&sm);
	elemOpAddTo(0, 1, 1.0, &sm);
	elemOpCreate(&sp);
	elemOpAddTo(1, 0, 1.0, &sp);
	elemOpCreate(&sz);
	elemOpAddTo(1, 1, 1.0, &sz);
	elemOpAddTo(0, 0, -1.0, &sz);

	hSingleAtom = prodSpaceCreate(2);
	hAtoms = prodSpaceCreate(0);
	for (int i = 0; i < nAtoms; ++i) {
		prodSpaceMul(hSingleAtom, &hAtoms);
	}
	hField = prodSpaceCreate(nPhotons + 1);
	hTot = prodSpaceCreate(0);
	prodSpaceMul(hAtoms, &hTot);
	prodSpaceMul(hField, &hTot);

	tensorOpNull(hAtoms, &inhomogeneousJz);
	for (int i = 0; i < nAtoms; ++i) {
		tensorOpAddScaledTo(omega(i), sz, i, inhomogeneousJz);
	}

	tensorOpNull(hTot, &jMinus);

	elemOpDestroy(&sm);
	elemOpDestroy(&sp);
	elemOpDestroy(&sz);
	prodSpaceDestroy(&hSingleAtom);
	prodSpaceDestroy(&hAtoms);
	prodSpaceDestroy(&hField);
	prodSpaceDestroy(&hTot);
	tensorOpDestroy(&inhomogeneousJz);
	tensorOpDestroy(&jMinus);

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}

