#include <stdlib.h>
#include <string.h>

#include "MboProdSpace.h"

struct MboProdSpace {
	int numSpaces;
	MboLocInd *dims;
};

MboProdSpace mboProdSpaceCreate(MboLocInd d)
{
	MboProdSpace sp = malloc(sizeof(*sp));
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

void mboProdSpaceDestroy(MboProdSpace *sp)
{
	free((*sp)->dims);
	free(*sp);
	*sp = 0;
}

void mboProdSpaceMul(MboProdSpace a, MboProdSpace *b)
{
	if (a->numSpaces == 0) return;
	(*b)->dims = realloc((*b)->dims, ((*b)->numSpaces + a->numSpaces) *
					     sizeof(*(*b)->dims));
	memmove((*b)->dims + a->numSpaces, (*b)->dims,
		(*b)->numSpaces * sizeof(*(*b)->dims));
	memcpy((*b)->dims, a->dims, a->numSpaces * sizeof(*(*b)->dims));
	(*b)->numSpaces += a->numSpaces;
}

MboProdSpace mboProdSpaceCopy(MboProdSpace sp)
{
	MboProdSpace copy = malloc(sizeof(*copy));
	copy->numSpaces = sp->numSpaces;
	copy->dims = malloc(sp->numSpaces * sizeof(*copy->dims));
	memcpy(copy->dims, sp->dims, sp->numSpaces * sizeof(*copy->dims));
	return copy;
}

MboGlobInd mboProdSpaceDim(MboProdSpace sp)
{
	MboGlobInd dim = 1;
	int i;
	for (i = 0; i < sp->numSpaces; ++i) {
		dim *= (MboGlobInd)sp->dims[i];
	}
	return dim;
}

int mboProdSpaceSize(MboProdSpace h)
{
	return h->numSpaces;
}

void mboProdSpaceGetDims(MboProdSpace h, int n, MboLocInd *dims)
{
	int i = 0;
	while (i < h->numSpaces && i < n) {
		dims[i] = h->dims[i];
		++i;
	}
}

int mboProdSpaceEqual(MboProdSpace h1, MboProdSpace h2)
{
	int i;
	if (mboProdSpaceSize(h1) != mboProdSpaceSize(h2)) return 0;
	for (i = 0; i < mboProdSpaceSize(h1); ++i) {
		if (h1->dims[i] != h2->dims[i]) return 0;
	}
	return 1;
}

int mboProdSpaceCheck(MboProdSpace h)
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

static int testMboProdSpaceCreate()
{
	int errs = 0;
	MboProdSpace sp;
	sp = mboProdSpaceCreate(2);
	CHK_EQUAL(sp->numSpaces, 1, errs);
	CHK_TRUE(sp->dims != 0, errs);
	mboProdSpaceDestroy(&sp);

	sp = mboProdSpaceCreate(0);
	CHK_EQUAL(sp->numSpaces, 0, errs);
	CHK_EQUAL(sp->dims, 0, errs);
	mboProdSpaceDestroy(&sp);

	return errs;
}

static int testMboProdSpaceDestroy()
{
	int errs = 0;
	MboProdSpace sp;

	sp = mboProdSpaceCreate(1);
	mboProdSpaceDestroy(&sp);
	return errs;
}

static int testMboProdSpaceMul()
{
	int errs = 0;
	MboProdSpace sp1;
	MboProdSpace sp2;

	sp1 = mboProdSpaceCreate(0);
	sp2 = mboProdSpaceCreate(0);
	mboProdSpaceMul(sp1, &sp2);
	CHK_EQUAL(mboProdSpaceDim(sp2), 1, errs);
	CHK_EQUAL(sp2->numSpaces, 0, errs);
	CHK_EQUAL(sp2->dims, 0, errs);
	mboProdSpaceDestroy(&sp1);
	mboProdSpaceDestroy(&sp2);

	sp1 = mboProdSpaceCreate(3);
	sp2 = mboProdSpaceCreate(2);
	CHK_EQUAL(mboProdSpaceDim(sp1), 3, errs);
	CHK_EQUAL(mboProdSpaceDim(sp2), 2, errs);
	mboProdSpaceMul(sp1, &sp2);
	CHK_EQUAL(mboProdSpaceDim(sp1), 3, errs);
	CHK_EQUAL(mboProdSpaceDim(sp2), 6, errs);
	CHK_EQUAL(sp2->numSpaces, 2, errs);
	CHK_EQUAL(sp2->dims[0], 3, errs);
	CHK_EQUAL(sp2->dims[1], 2, errs);
	mboProdSpaceDestroy(&sp1);
	mboProdSpaceDestroy(&sp2);

	sp1 = mboProdSpaceCreate(0);
	sp2 = mboProdSpaceCreate(2);
	mboProdSpaceMul(sp1, &sp2);
	CHK_EQUAL(mboProdSpaceDim(sp2), 2, errs);
	CHK_EQUAL(sp2->numSpaces, 1, errs);
	CHK_EQUAL(sp2->dims[0], 2, errs);
	mboProdSpaceDestroy(&sp1);
	mboProdSpaceDestroy(&sp2);

	sp1 = mboProdSpaceCreate(5);
	sp2 = mboProdSpaceCreate(0);
	mboProdSpaceMul(sp1, &sp2);
	CHK_EQUAL(mboProdSpaceDim(sp2), 5, errs);
	CHK_EQUAL(sp2->numSpaces, 1, errs);
	CHK_EQUAL(sp2->dims[0], 5, errs);
	mboProdSpaceDestroy(&sp1);
	mboProdSpaceDestroy(&sp2);

	return errs;
}

static int testBuildSpace()
{
	int errs = 0;
	int N = 20;
	int i;
	MboProdSpace h;
	MboProdSpace hTot;

	h = mboProdSpaceCreate(2);
	CHK_EQUAL(mboProdSpaceDim(h), 2, errs);
	hTot = mboProdSpaceCreate(0);
	CHK_EQUAL(mboProdSpaceDim(hTot), 1, errs);
	for (i = 0; i < N; ++i) {
		mboProdSpaceMul(h, &hTot);
	}
	CHK_EQUAL(mboProdSpaceDim(hTot), (MboGlobInd)1 << N, errs);

	mboProdSpaceDestroy(&h);
	mboProdSpaceDestroy(&hTot);
	return errs;
}

static int testMultiplyWithSelf()
{
	int errs = 0;
	int i;
        long long d;
	MboProdSpace h;

	h = mboProdSpaceCreate(2);
        d = mboProdSpaceDim(h);
	CHK_EQUAL(mboProdSpaceDim(h), 2, errs);
	for (i = 0; i < 3; ++i) {
		d *= d;
		mboProdSpaceMul(h, &h);
		CHK_EQUAL(mboProdSpaceDim(h), d, errs);
	}
        mboProdSpaceDestroy(&h);
	return errs;
}

static int testMultiplyLargeDims()
{
	int errs = 0;
        MboLocInd d1, d2;
	MboProdSpace h1;
	MboProdSpace h2;

        d1 = 0x7FFFFFFF;
        d2 = 0x7FFFFFFF;
	h1 = mboProdSpaceCreate(d1);
	h2 = mboProdSpaceCreate(d2);
	mboProdSpaceMul(h1, &h2);
	CHK_EQUAL(mboProdSpaceDim(h2), (MboGlobInd)d1 * (MboGlobInd)d2, errs);
	mboProdSpaceDestroy(&h1);
        mboProdSpaceDestroy(&h2);

	return errs;
}

static int testMboProdSpaceSize()
{
	MboProdSpace h, h2;
	int s, errs = 0;

	h = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCreate(6);
	s = mboProdSpaceSize(h);
	CHK_EQUAL(s, 1, errs);

	mboProdSpaceMul(h2, &h);
	mboProdSpaceMul(h2, &h);
	s = mboProdSpaceSize(h);
	CHK_EQUAL(s, 3, errs);

	mboProdSpaceDestroy(&h);
	mboProdSpaceDestroy(&h2);
	return errs;
}

static int testMboProdSpaceGetDims()
{
	int errs = 0;
	MboLocInd dims[10];
	MboProdSpace h;

	h = mboProdSpaceCreate(0);
	dims[0] = 7;
	mboProdSpaceGetDims(h, 0, dims);
	CHK_EQUAL(dims[0], 7, errs);
	mboProdSpaceGetDims(h, 2, dims);
	CHK_EQUAL(dims[0], 7, errs);
	mboProdSpaceDestroy(&h);

	h = mboProdSpaceCreate(20);
	dims[0] = 7;
	dims[1] = 13;
	mboProdSpaceGetDims(h, 0, dims);
	CHK_EQUAL(dims[0], 7, errs);
	mboProdSpaceGetDims(h, 2, dims);
	CHK_EQUAL(dims[0], 20, errs);
	CHK_EQUAL(dims[1], 13, errs);
	mboProdSpaceDestroy(&h);

	return errs;
}

static int testMboProdSpaceEqual()
{
	MboProdSpace h1, h2, h3;
	int errs = 0;

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCopy(h1);
	CHK_TRUE(mboProdSpaceEqual(h1, h2), errs);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCopy(h1);
	mboProdSpaceMul(h1, &h2);
	CHK_FALSE(mboProdSpaceEqual(h1, h2), errs);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCopy(h1);
	mboProdSpaceMul(h1, &h2);
	mboProdSpaceMul(h2, &h1);
	CHK_FALSE(mboProdSpaceEqual(h1, h2), errs);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);

	h1 = mboProdSpaceCreate(2);
	h2 = mboProdSpaceCopy(h1);
	h3 = mboProdSpaceCopy(h2);
	mboProdSpaceMul(h1, &h2);
	mboProdSpaceMul(h1, &h2);
	mboProdSpaceMul(h1, &h3);
	mboProdSpaceMul(h1, &h3);
	CHK_TRUE(mboProdSpaceEqual(h2, h3), errs);
	mboProdSpaceDestroy(&h1);
	mboProdSpaceDestroy(&h2);
	mboProdSpaceDestroy(&h3);

	return errs;
}

static int testMboProdSpaceCheck()
{
	MboProdSpace h, h2;
	int errs = 0;

	h = mboProdSpaceCreate(2);
	mboProdSpaceCheck(h);
	h2 = mboProdSpaceCreate(0);
	mboProdSpaceCheck(h2);

	mboProdSpaceMul(h, &h2);
	mboProdSpaceMul(h, &h2);
	mboProdSpaceCheck(h2);

	mboProdSpaceDestroy(&h);
	mboProdSpaceDestroy(&h2);
	return errs;
}

int mboProdSpaceTest()
{
	int errs = 0;
	errs += testMboProdSpaceCreate();
	errs += testMboProdSpaceDestroy();
	errs += testMboProdSpaceMul();
	errs += testBuildSpace();
	errs += testMultiplyWithSelf();
	errs += testMultiplyLargeDims();
	errs += testMboProdSpaceSize();
	errs += testMboProdSpaceGetDims();
	errs += testMboProdSpaceEqual();
	errs += testMboProdSpaceCheck();
	return errs;
}
