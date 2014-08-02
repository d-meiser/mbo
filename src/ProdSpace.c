#include <stdlib.h>
#include <string.h>

#include "ProdSpace.h"

struct ProdSpace {
	int numSpaces;
	int *dims;
};

ProdSpace prodSpaceCreate(int d)
{
	ProdSpace sp = malloc(sizeof(*sp));
	if (d) {
		sp->numSpaces = 1;
		sp->dims = malloc(1 * sizeof(*sp->dims));
		sp->dims[0] = d;
	} else {
		sp->numSpaces = 0;
		sp->dims = 0;
	}
	return sp;
}

void prodSpaceDestroy(ProdSpace *sp)
{
	free((*sp)->dims);
	free(*sp);
	*sp = 0;
}

void prodSpaceMul(ProdSpace a, ProdSpace *b)
{
	if (a->numSpaces == 0) return;
	(*b)->dims = realloc((*b)->dims, ((*b)->numSpaces + a->numSpaces) *
					     sizeof(*(*b)->dims));
	memmove((*b)->dims + a->numSpaces, (*b)->dims,
		(*b)->numSpaces * sizeof(*(*b)->dims));
	memcpy((*b)->dims, a->dims, a->numSpaces * sizeof(*(*b)->dims));
	(*b)->numSpaces += a->numSpaces;
}

ProdSpace prodSpaceCopy(ProdSpace sp)
{
	ProdSpace copy = malloc(sizeof(*copy));
	copy->numSpaces = sp->numSpaces;
	copy->dims = malloc(sp->numSpaces * sizeof(*copy->dims));
	memcpy(copy->dims, sp->dims, sp->numSpaces * sizeof(*copy->dims));
	return copy;
}

long long prodSpaceDim(ProdSpace sp)
{
	long long dim = 1;
	int i;
	for (i = 0; i < sp->numSpaces; ++i) {
		dim *= (long long)sp->dims[i];
	}
	return dim;
}

int prodSpaceSize(ProdSpace h)
{
	return h->numSpaces;
}

void prodSpaceGetDims(ProdSpace h, int n, int *dims)
{
	int i = 0;
	while (i < h->numSpaces && i < n) {
		dims[i] = h->dims[i];
		++i;
	}
}

int prodSpaceEqual(ProdSpace h1, ProdSpace h2)
{
	int i;
	if (prodSpaceSize(h1) != prodSpaceSize(h2)) return 0;
	for (i = 0; i < prodSpaceSize(h1); ++i) {
		if (h1->dims[i] != h2->dims[i]) return 0;
	}
	return 1;
}

int prodSpaceCheck(ProdSpace h)
{
	int errs = 0, i;
	for (i = 0; i < h->numSpaces; ++i) {
		if (h->dims[i] <= 0) ++errs;
	}
	return errs;
}

/*
 * Tests
 * */
#include "TestUtils.h"

static int testprodSpaceCreate()
{
	int errs = 0;
	ProdSpace sp;
	sp = prodSpaceCreate(2);
	CHK_EQUAL(sp->numSpaces, 1, errs);
	CHK_TRUE(sp->dims != 0, errs);
	prodSpaceDestroy(&sp);

	sp = prodSpaceCreate(0);
	CHK_EQUAL(sp->numSpaces, 0, errs);
	CHK_EQUAL(sp->dims, 0, errs);
	prodSpaceDestroy(&sp);

	return errs;
}

static int testprodSpaceDestroy()
{
	int errs = 0;
	ProdSpace sp;

	sp = prodSpaceCreate(1);
	prodSpaceDestroy(&sp);
	return errs;
}

static int testprodSpaceMul()
{
	int errs = 0;
	ProdSpace sp1;
	ProdSpace sp2;

	sp1 = prodSpaceCreate(0);
	sp2 = prodSpaceCreate(0);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(prodSpaceDim(sp2), 1, errs);
	CHK_EQUAL(sp2->numSpaces, 0, errs);
	CHK_EQUAL(sp2->dims, 0, errs);
	prodSpaceDestroy(&sp1);
	prodSpaceDestroy(&sp2);

	sp1 = prodSpaceCreate(3);
	sp2 = prodSpaceCreate(2);
	CHK_EQUAL(prodSpaceDim(sp1), 3, errs);
	CHK_EQUAL(prodSpaceDim(sp2), 2, errs);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(prodSpaceDim(sp1), 3, errs);
	CHK_EQUAL(prodSpaceDim(sp2), 6, errs);
	CHK_EQUAL(sp2->numSpaces, 2, errs);
	CHK_EQUAL(sp2->dims[0], 3, errs);
	CHK_EQUAL(sp2->dims[1], 2, errs);
	prodSpaceDestroy(&sp1);
	prodSpaceDestroy(&sp2);

	sp1 = prodSpaceCreate(0);
	sp2 = prodSpaceCreate(2);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(prodSpaceDim(sp2), 2, errs);
	CHK_EQUAL(sp2->numSpaces, 1, errs);
	CHK_EQUAL(sp2->dims[0], 2, errs);
	prodSpaceDestroy(&sp1);
	prodSpaceDestroy(&sp2);

	sp1 = prodSpaceCreate(5);
	sp2 = prodSpaceCreate(0);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(prodSpaceDim(sp2), 5, errs);
	CHK_EQUAL(sp2->numSpaces, 1, errs);
	CHK_EQUAL(sp2->dims[0], 5, errs);
	prodSpaceDestroy(&sp1);
	prodSpaceDestroy(&sp2);

	return errs;
}

static int testBuildSpace()
{
	int errs = 0;
	int N = 20;
	int i;
	ProdSpace h;
	ProdSpace hTot;

	h = prodSpaceCreate(2);
	CHK_EQUAL(prodSpaceDim(h), 2, errs);
	hTot = prodSpaceCreate(0);
	CHK_EQUAL(prodSpaceDim(hTot), 1, errs);
	for (i = 0; i < N; ++i) {
		prodSpaceMul(h, &hTot);
	}
	CHK_EQUAL(prodSpaceDim(hTot), 1 << N, errs);

	prodSpaceDestroy(&h);
	prodSpaceDestroy(&hTot);
	return errs;
}

static int testMultiplyWithSelf()
{
	int errs = 0;
	int i;
        long long d;
	ProdSpace h;

	h = prodSpaceCreate(2);
        d = prodSpaceDim(h);
	CHK_EQUAL(prodSpaceDim(h), 2, errs);
	for (i = 0; i < 3; ++i) {
		d *= d;
		prodSpaceMul(h, &h);
		CHK_EQUAL(prodSpaceDim(h), d, errs);
	}
        prodSpaceDestroy(&h);
	return errs;
}

static int testMultiplyLargeDims()
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
	CHK_EQUAL(prodSpaceDim(h2), d1 * d2, errs);
	prodSpaceDestroy(&h1);
        prodSpaceDestroy(&h2);

	return errs;
}

static int testSize()
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

static int testGetDims()
{
	int errs = 0;
	int dims[10];
	ProdSpace h;

	h = prodSpaceCreate(0);
	dims[0] = 7;
	prodSpaceGetDims(h, 0, dims);
	CHK_EQUAL(dims[0], 7, errs);
	prodSpaceGetDims(h, 2, dims);
	CHK_EQUAL(dims[0], 7, errs);
	prodSpaceDestroy(&h);

	h = prodSpaceCreate(20);
	dims[0] = 7;
	dims[1] = 13;
	prodSpaceGetDims(h, 0, dims);
	CHK_EQUAL(dims[0], 7, errs);
	prodSpaceGetDims(h, 2, dims);
	CHK_EQUAL(dims[0], 20, errs);
	CHK_EQUAL(dims[1], 13, errs);
	prodSpaceDestroy(&h);

	return errs;
}

static int testEqual()
{
	ProdSpace h1, h2, h3;
	int errs = 0;

	h1 = prodSpaceCreate(2);
	h2 = prodSpaceCopy(h1);
	CHK_TRUE(prodSpaceEqual(h1, h2), errs);
	prodSpaceDestroy(&h1);
	prodSpaceDestroy(&h2);

	h1 = prodSpaceCreate(2);
	h2 = prodSpaceCopy(h1);
	prodSpaceMul(h1, &h2);
	CHK_FALSE(prodSpaceEqual(h1, h2), errs);
	prodSpaceDestroy(&h1);
	prodSpaceDestroy(&h2);

	h1 = prodSpaceCreate(2);
	h2 = prodSpaceCopy(h1);
	prodSpaceMul(h1, &h2);
	prodSpaceMul(h2, &h1);
	CHK_FALSE(prodSpaceEqual(h1, h2), errs);
	prodSpaceDestroy(&h1);
	prodSpaceDestroy(&h2);

	h1 = prodSpaceCreate(2);
	h2 = prodSpaceCopy(h1);
	h3 = prodSpaceCopy(h2);
	prodSpaceMul(h1, &h2);
	prodSpaceMul(h1, &h2);
	prodSpaceMul(h1, &h3);
	prodSpaceMul(h1, &h3);
	CHK_TRUE(prodSpaceEqual(h2, h3), errs);
	prodSpaceDestroy(&h1);
	prodSpaceDestroy(&h2);
	prodSpaceDestroy(&h3);

	return errs;
}

static int testCheck()
{
	ProdSpace h, h2;
	int errs = 0;

	h = prodSpaceCreate(2);
	prodSpaceCheck(h);
	h2 = prodSpaceCreate(0);
	prodSpaceCheck(h2);

	prodSpaceMul(h, &h2);
	prodSpaceMul(h, &h2);
	prodSpaceCheck(h2);

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
	errs += testGetDims();
	errs += testEqual();
	errs += testCheck();
	return errs;
}
