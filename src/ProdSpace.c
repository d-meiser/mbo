#include <stdlib.h>

#include "ProdSpace.h"

struct ProdSpace {
	int dim;
	struct ProdSpace *next;
};

ProdSpace prodSpaceCreate(int d)
{
	ProdSpace sp = malloc(sizeof(*sp));
	sp->dim = d;
	sp->next = 0;
	return sp;
}

void prodSpaceDestroy(ProdSpace *sp)
{
	if (*sp == 0) return;
	if ((*sp)->next) {
		prodSpaceDestroy(&(*sp)->next);
	}
	free(*sp);
	*sp = 0;
}

void prodSpaceMul(ProdSpace a, ProdSpace *b)
{
	if (a->dim == 0) return;
	if ((*b)->dim == 0) {
		free(*b);
		*b = 0;
	}
	ProdSpace aCopy = prodSpaceCopy(a);
	ProdSpace last = aCopy;
	while (last->next) {
		last = last->next;
	}
	last->next = *b;
	*b = aCopy;
}

ProdSpace prodSpaceCopy(ProdSpace sp)
{
	ProdSpace copy = 0;
	if (sp) {
		copy = malloc(sizeof(*copy));
		copy->dim = sp->dim;
		copy->next = prodSpaceCopy(sp->next);
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

int prodSpaceSize(ProdSpace h)
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

int testprodSpaceCreate()
{
	int errs = 0;
	ProdSpace sp;
	sp = prodSpaceCreate(2);
	prodSpaceDestroy(&sp);
	return errs;
}

int testprodSpaceDestroy()
{
	int errs = 0;
	ProdSpace sp;

	sp = prodSpaceCreate(0);
	prodSpaceDestroy(&sp);
	return errs;
}

int testprodSpaceMul()
{
	int errs = 0;
	ProdSpace sp1;
	ProdSpace sp2;

	sp1 = prodSpaceCreate(2);
	sp2 = prodSpaceCreate(2);

	CHK_EQUAL(DimProdSpace(sp1), 2, errs);
	CHK_EQUAL(DimProdSpace(sp2), 2, errs);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(DimProdSpace(sp1), 2, errs);
	CHK_EQUAL(DimProdSpace(sp2), 4, errs);

	prodSpaceDestroy(&sp1);
	prodSpaceDestroy(&sp2);
	return errs;
}

int testBuildSpace()
{
	int errs = 0;
	int N = 20;
	int i;
	ProdSpace h;
	ProdSpace hTot;

	h = prodSpaceCreate(2);
	CHK_EQUAL(DimProdSpace(h), 2, errs);
	hTot = prodSpaceCreate(0);
	CHK_EQUAL(DimProdSpace(hTot), 0, errs);
	for (i = 0; i < N; ++i) {
		prodSpaceMul(h, &hTot);
	}
	CHK_EQUAL(DimProdSpace(hTot), 1 << N, errs);

	prodSpaceDestroy(&h);
	prodSpaceDestroy(&hTot);
	return errs;
}

int testMultiplyWithSelf()
{
	int errs = 0;
	int i;
        long long d;
	ProdSpace h;

	h = prodSpaceCreate(2);
        d = DimProdSpace(h);
	CHK_EQUAL(DimProdSpace(h), 2, errs);
	for (i = 0; i < 3; ++i) {
		d *= d;
		prodSpaceMul(h, &h);
		CHK_EQUAL(DimProdSpace(h), d, errs);
	}
        prodSpaceDestroy(&h);
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
	h1 = prodSpaceCreate(d1);
	h2 = prodSpaceCreate(d2);
	prodSpaceMul(h1, &h2);
	CHK_EQUAL(DimProdSpace(h2), d1 * d2, errs);
	prodSpaceDestroy(&h1);
        prodSpaceDestroy(&h2);

	return errs;
}

int testSize()
{
	ProdSpace h, h2;
	int s, errs = 0;

	h = prodSpaceCreate(2);
	h2 = prodSpaceCreate(6);
	s = prodSpaceSize(h);
	CHK_EQUAL(s, 1, errs);

	prodSpaceMul(h2, &h);
	prodSpaceMul(h2, &h);
	s = prodSpaceSize(h);
	CHK_EQUAL(s, 3, errs);

	prodSpaceDestroy(&h);
	prodSpaceDestroy(&h2);
	return errs;
}

int prodSpaceTest()
{
	int errs = 0;
	errs += testprodSpaceCreate();
	errs += testprodSpaceDestroy();
	errs += testprodSpaceMul();
	errs += testBuildSpace();
        errs += testMultiplyWithSelf();
        errs += testMultiplyLargeDims();
	errs += testSize();
	return errs;
}
