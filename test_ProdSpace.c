#include "ProdSpace.h"

int testCreateProdSpace()
{
	int errs = 0;
	ProdSpace sp;
	sp = CreateProdSpace(2);
	DestroyProdSpace(sp);
	return errs;
}

int testMultToProdSpace()
{
	int errs = 0;
	ProdSpace sp1;
	ProdSpace sp2;

	sp1 = CreateProdSpace(2);
	sp2 = CreateProdSpace(2);

	MultToProdSpace(sp1, sp2);

	DestroyProdSpace(sp1);
	DestroyProdSpace(sp2);
	return errs;
}

int main()
{
	int errs = 0;
	errs += testCreateProdSpace();
	errs += testMultToProdSpace();
	return errs;
}
