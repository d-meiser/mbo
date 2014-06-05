#include <Quo.h>
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
	TOp inhomogeneousJz;
	TOp jMinus;

	CreateElemOp(&sm);
	AddToElemOp(0, 1, 1.0, &sm);
	CreateElemOp(&sp);
	AddToElemOp(1, 0, 1.0, &sp);
	CreateElemOp(&sz);
	AddToElemOp(1, 1, 1.0, &sz);
	AddToElemOp(0, 0, -1.0, &sz);

	hSingleAtom = CreateProdSpace(2);
	hAtoms = CreateProdSpace(0);
	for (int i = 0; i < nAtoms; ++i) {
		MultToProdSpace(hSingleAtom, &hAtoms);
	}
	hField = CreateProdSpace(nPhotons + 1);
	hTot = CreateProdSpace(0);
	MultToProdSpace(hAtoms, &hTot);
	MultToProdSpace(hField, &hTot);

	CreateTOp(hAtoms, &inhomogeneousJz);
	for (int i = 0; i < nAtoms; ++i) {
		AddScaledToTOp(omega(i), sz, i, inhomogeneousJz);
	}

	CreateTOp(hTot, &jMinus);

	DestroyElemOp(&sm);
	DestroyElemOp(&sp);
	DestroyElemOp(&sz);
	DestroyProdSpace(&hSingleAtom);
	DestroyProdSpace(&hAtoms);
	DestroyProdSpace(&hField);
	DestroyProdSpace(&hTot);
	DestroyTOp(&inhomogeneousJz);
	DestroyTOp(&jMinus);

	return 0;
}

double omega(int i)
{
	return 1.3 * i;
}

