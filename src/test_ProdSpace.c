#include "ProdSpace.h"
#include <stdio.h>

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

	MultToProdSpace(sp1, &sp2);

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
        errs += (2 != DimProdSpace(h));
        printf("DimProdSpace(h) == %lld\n", DimProdSpace(h));
        hTot = CreateProdSpace(0);
        errs += (0 != DimProdSpace(hTot));
        printf("DimProdSpace(hTot) == %lld\n", DimProdSpace(hTot));
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &hTot);
	}
        errs += ((1 << N) != DimProdSpace(hTot));
        printf("DimProdSpace(hTot) == %lld\n", DimProdSpace(hTot));

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
