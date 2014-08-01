#include <stdlib.h>
#include <string.h>

#include "ProdSpace.h"

struct ProdSpace {
	int numSpaces;
	int numRefs;
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
	sp->numRefs = 1;
	return sp;
}

void prodSpaceRetain(ProdSpace sp)
{
	++sp->numRefs;
}

void prodSpaceRelease(ProdSpace sp)
{
	--sp->numRefs;
	if (!sp->numRefs) {
		free(sp->dims);
		free(sp);
	}
}

void prodSpaceMul(ProdSpace a, ProdSpace *b)
{
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

int testprodSpaceCreate()
{
	int errs = 0;
	ProdSpace sp;
	sp = prodSpaceCreate(2);
	prodSpaceRelease(sp);
	return errs;
}

int testprodSpaceRetain()
{
	int errs = 0;
	ProdSpace sp;

	sp = prodSpaceCreate(1);
	prodSpaceRetain(sp);
	prodSpaceRetain(sp);
	prodSpaceRelease(sp);
	prodSpaceRelease(sp);
	prodSpaceRelease(sp);
	return errs;
}

int testprodSpaceRelease()
{
	int errs = 0;
	ProdSpace sp;

	sp = prodSpaceCreate(1);
	prodSpaceRelease(sp);
	return errs;
}

int testprodSpaceMul()
{
	int errs = 0;
	ProdSpace sp1;
	ProdSpace sp2;

	sp1 = prodSpaceCreate(2);
	sp2 = prodSpaceCreate(2);

	CHK_EQUAL(prodSpaceDim(sp1), 2, errs);
	CHK_EQUAL(prodSpaceDim(sp2), 2, errs);
	prodSpaceMul(sp1, &sp2);
	CHK_EQUAL(prodSpaceDim(sp1), 2, errs);
	CHK_EQUAL(prodSpaceDim(sp2), 4, errs);

	prodSpaceRelease(sp1);
	prodSpaceRelease(sp2);
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
	CHK_EQUAL(prodSpaceDim(h), 2, errs);
	hTot = prodSpaceCreate(0);
	CHK_EQUAL(prodSpaceDim(hTot), 1, errs);
	for (i = 0; i < N; ++i) {
		prodSpaceMul(h, &hTot);
	}
	CHK_EQUAL(prodSpaceDim(hTot), 1 << N, errs);

	prodSpaceRelease(h);
	prodSpaceRelease(hTot);
	return errs;
}

int testMultiplyWithSelf()
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
        prodSpaceRelease(h);
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
	CHK_EQUAL(prodSpaceDim(h2), d1 * d2, errs);
	prodSpaceRelease(h1);
        prodSpaceRelease(h2);

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

	prodSpaceRelease(h);
	prodSpaceRelease(h2);
	return errs;
}

int prodSpaceTest()
{
	int errs = 0;
	errs += testprodSpaceCreate();
	errs += testprodSpaceRetain();
	errs += testprodSpaceRelease();
	errs += testprodSpaceMul();
	errs += testBuildSpace();
        errs += testMultiplyWithSelf();
        errs += testMultiplyLargeDims();
	errs += testSize();
	return errs;
}
