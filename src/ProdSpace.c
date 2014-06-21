#include <stdlib.h>

#include "ProdSpace.h"

struct ProdSpace {
	int dim;
	struct ProdSpace *next;
};

ProdSpace CreateProdSpace(int d)
{
	ProdSpace sp = malloc(sizeof(*sp));
	sp->dim = d;
	sp->next = 0;
	return sp;
}

void DestroyProdSpace(ProdSpace *sp)
{
	if (*sp == 0) return;
	if ((*sp)->next) {
		DestroyProdSpace(&(*sp)->next);
	}
	free(*sp);
	*sp = 0;
}

void MultToProdSpace(ProdSpace a, ProdSpace *b)
{
	if (a->dim == 0) return;
	if ((*b)->dim == 0) {
		free(*b);
		*b = 0;
	}
	ProdSpace aCopy = CopyProdSpace(a);
	ProdSpace last = aCopy;
	while (last->next) {
		last = last->next;
	}
	last->next = *b;
	*b = aCopy;
}

ProdSpace CopyProdSpace(ProdSpace sp)
{
	ProdSpace copy = 0;
	if (sp) {
		copy = malloc(sizeof(*copy));
		copy->dim = sp->dim;
		copy->next = CopyProdSpace(sp->next);
	}
	return copy;
}

long long DimProdSpace(ProdSpace sp)
{
	long long dim = 1;
	while (sp) {
		dim *= sp->dim;
		sp = sp->next;
	}
	return dim;
}

int SizeProdSpace(ProdSpace h)
{
	int size = 0;
	for (; h != 0; h = h->next) {
		++size;
	}
	return size;
}

int prodSpaceCheck(ProdSpace h)
{
	int errs = 0;
	while (h) {
		if (h->dim < 0) ++errs;
		h = h->next;
	}
	return errs;
}

/* 
 * Tests
 * */
#include "TestUtils.h"

int testCreateProdSpace()
{
	int errs = 0;
	ProdSpace sp;
	sp = CreateProdSpace(2);
	DestroyProdSpace(&sp);
	return errs;
}

int testDestroyProdSpace()
{
	int errs = 0;
	ProdSpace sp;

	sp = CreateProdSpace(0);
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

int testMultiplyWithSelf()
{
	int errs = 0;
	int i;
        long long d;
	ProdSpace h;

	h = CreateProdSpace(2);
        d = DimProdSpace(h);
	CHK_EQUAL(DimProdSpace(h), 2, errs);
	for (i = 0; i < 3; ++i) {
		d *= d;
		MultToProdSpace(h, &h);
		CHK_EQUAL(DimProdSpace(h), d, errs);
	}
        DestroyProdSpace(&h);
	return errs;
}

int testMultiplyLargeDims()
{
	int errs = 0;
        long long d1, d2;
	ProdSpace h1;
	ProdSpace h2;

        d1 = 0x7FFFFFFF;
        d2 = 0x7FFFFFFF;
	h1 = CreateProdSpace(d1);
	h2 = CreateProdSpace(d2);
	MultToProdSpace(h1, &h2);
	CHK_EQUAL(DimProdSpace(h2), d1 * d2, errs);
	DestroyProdSpace(&h1);
        DestroyProdSpace(&h2);

	return errs;
}

int testSize()
{
	ProdSpace h, h2;
	int s, errs = 0;

	h = CreateProdSpace(2);
	h2 = CreateProdSpace(6);
	s = SizeProdSpace(h);
	CHK_EQUAL(s, 1, errs);

	MultToProdSpace(h2, &h);
	MultToProdSpace(h2, &h);
	s = SizeProdSpace(h);
	CHK_EQUAL(s, 3, errs);

	DestroyProdSpace(&h);
	DestroyProdSpace(&h2);
	return errs;
}

int prodSpaceTest()
{
	int errs = 0;
	errs += testCreateProdSpace();
	errs += testDestroyProdSpace();
	errs += testMultToProdSpace();
	errs += testBuildSpace();
        errs += testMultiplyWithSelf();
        errs += testMultiplyLargeDims();
	errs += testSize();
	return errs;
}
