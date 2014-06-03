#include "ProdSpace.h"
#include <stdio.h>

#include "TestUtils.h"

int testCreateProdSpace()
{
	int errs = 0;
	ProdSpace sp;
	sp = CreateProdSpace(2);
	DestroyProdSpace(&sp);
	return errs;
}

int testMultToProdSpace()
{
	int errs = 0;
	ProdSpace sp1;
	ProdSpace sp2;

	sp1 = CreateProdSpace(2);
	sp2 = CreateProdSpace(2);

	CHK_EQUAL(DimProdSpace(sp1), 2, errs);
	CHK_EQUAL(DimProdSpace(sp2), 2, errs);
	MultToProdSpace(sp1, &sp2);
	CHK_EQUAL(DimProdSpace(sp1), 2, errs);
	CHK_EQUAL(DimProdSpace(sp2), 4, errs);

	DestroyProdSpace(&sp1);
	DestroyProdSpace(&sp2);
	return errs;
}

int testBuildSpace()
{
	int errs = 0;
	int N = 20;
	int i;
	ProdSpace h;
	ProdSpace hTot;

	h = CreateProdSpace(2);
	CHK_EQUAL(DimProdSpace(h), 2, errs);
	hTot = CreateProdSpace(0);
	CHK_EQUAL(DimProdSpace(hTot), 0, errs);
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &hTot);
	}
	CHK_EQUAL(DimProdSpace(hTot), 1 << N, errs);

	DestroyProdSpace(&h);
	DestroyProdSpace(&hTot);
	return errs;
}

int main()
{
	int errs = 0;
	errs += testCreateProdSpace();
	errs += testMultToProdSpace();
	errs += testBuildSpace();
	return errs;
}
