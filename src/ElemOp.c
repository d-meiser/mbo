#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ElemOp.h"
#include "NonZeroEntry.h"

struct ElemOp
{
	int nOps;
	struct NonZeroEntry *entries;
};

void elemOpCreate(ElemOp *eo)
{
	*eo = malloc(sizeof(struct ElemOp));
	(*eo)->nOps = 0;
	(*eo)->entries = 0;
}

void elemOpDestroy(ElemOp *eo)
{
	free((*eo)->entries);
	free(*eo);
	*eo = 0;
}

void elemOpAddTo(int m, int n, double val, ElemOp *eo)
{
	ElemOp a = *eo;
	a->entries = realloc(a->entries, (a->nOps + 1) * sizeof(*a->entries));
	a->entries[a->nOps].m = m;
	a->entries[a->nOps].n = n;
	a->entries[a->nOps].val = val;
	++a->nOps;
}

void elemOpScale(double alpha, ElemOp eo)
{
	int i;
	for (i = 0; i < eo->nOps; ++i) {
		eo->entries[i].val *= alpha;
	}
}

void elemOpPlus(ElemOp a, ElemOp *b)
{
	(*b)->entries = realloc((*b)->entries, ((*b)->nOps + a->nOps) *
						   sizeof(*(*b)->entries));
	memcpy((*b)->entries + (*b)->nOps, a->entries,
	       a->nOps * sizeof(*a->entries));
	(*b)->nOps += a->nOps;
}

int elemOpCheck(ElemOp a)
{
	int errs = 0;
	int i;
	for (i = 0; i < a->nOps; ++i) {
		if (a->entries[i].m < 0) ++errs;
		if (a->entries[i].n < 0) ++errs;
	}
	return errs;
}

ElemOp sigmaPlus()
{
	ElemOp sp;
	elemOpCreate(&sp);
	elemOpAddTo(1, 0, 1.0, &sp);
	return sp;
}

ElemOp sigmaMinus()
{
	ElemOp sm;
	elemOpCreate(&sm);
	elemOpAddTo(0, 1, 1.0, &sm);
	return sm;
}

ElemOp sigmaZ()
{
	ElemOp sz;
	elemOpCreate(&sz);
	elemOpAddTo(1, 1, 1.0, &sz);
	elemOpAddTo(0, 0, -1.0, &sz);
	return sz;
}

ElemOp eye(int d)
{
	ElemOp e;
	int i;
	elemOpCreate(&e);
	for (i = 0; i < d; ++i) {
		elemOpAddTo(i, i, 1.0, &e);
	}
	return e;
}

ElemOp numOp(int d)
{
	ElemOp n;
	int i;
	elemOpCreate(&n);
	for (i = 1; i < d; ++i) {
		elemOpAddTo(i, i, i, &n);
	}
	return n;
}

ElemOp annihilationOp(int d)
{
	ElemOp a;
	int i;
	elemOpCreate(&a);
	for (i = 1; i < d; ++i) {
		elemOpAddTo(i - 1, i, sqrt(i), &a);
	}
	return a;
}

ElemOp creationOp(int d)
{
	ElemOp ad;
	int i;
	elemOpCreate(&ad);
	for (i = 1; i < d; ++i) {
		elemOpAddTo(i, i - 1, sqrt(i), &ad);
	}
	return ad;
}

void elemOpMul(ElemOp a, ElemOp *b)
{
	ElemOp prod;
	int numOps, i, j;
	elemOpCreate(&prod);
	numOps = 0;
	for (i = 0; i < a->nOps; ++i) {
		for (j = 0; j < (*b)->nOps; ++j) {
			if (a->entries[i].n == (*b)->entries[j].m) {
				++numOps;
			}
		}
	}
	prod->entries = realloc(prod->entries, numOps * sizeof(*prod->entries));
	prod->nOps = numOps;
	numOps = 0;
	for (i = 0; i < a->nOps; ++i) {
		for (j = 0; j < (*b)->nOps; ++j) {
			if (a->entries[i].n == (*b)->entries[j].m) {
				prod->entries[numOps].m = a->entries[i].m;
				prod->entries[numOps].n = (*b)->entries[j].n;
				prod->entries[numOps].val =
				    a->entries[i].val * (*b)->entries[j].val;
				++numOps;
			}
		}
	}
	elemOpDestroy(b);
	*b = prod;
}

ElemOp elemOpCopy(ElemOp a)
{
	ElemOp copy;
	elemOpCreate(&copy);
	copy->nOps = a->nOps;
	copy->entries =
	    realloc(copy->entries, a->nOps * sizeof(*copy->entries));
	memcpy(copy->entries, a->entries, a->nOps * sizeof(*copy->entries));
	return copy;
}

/*
 * Tests
 * */

#include "TestUtils.h"
#define EPS 1.0e-12

static int testElemOpCreate()
{
	int errs = 0;
	ElemOp a;
	elemOpCreate(&a);

	CHK_EQUAL(0, a->nOps, errs);
	CHK_EQUAL(0, a->entries, errs);

	elemOpDestroy(&a);
	return errs;
}

static int testElemOpAddTo()
{
	int errs = 0;
	ElemOp a;
	elemOpCreate(&a);

	elemOpAddTo(1, 2, 3.0, &a);

	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 1, errs);
	CHK_EQUAL(a->entries[0].n, 2, errs);
	CHK_EQUAL(a->entries[0].val, 3.0, errs);

	elemOpDestroy(&a);
	return errs;
}

static int testElemOpScale()
{
	int errs = 0;
	ElemOp a;
	elemOpCreate(&a);

	elemOpAddTo(1, 2, 3.0, &a);
	elemOpScale(3.4, a);
	CHK_CLOSE(a->entries[0].val, 3.0 * 3.4, EPS, errs);

	elemOpDestroy(&a);
	return errs;
}

static int testElemOpPlus()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val, 3.0, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val, 3.0, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpAddTo(0, 1, 3.0, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->entries[1].m, 0, errs);
	CHK_EQUAL(b->entries[1].n, 1, errs);
	CHK_CLOSE(b->entries[1].val, 3.0, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpAddTo(5, 11, -2.0, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 5, errs);
	CHK_EQUAL(b->entries[0].n, 11, errs);
	CHK_CLOSE(b->entries[0].val, -2.0, EPS, errs);
	CHK_EQUAL(b->entries[1].m, 0, errs);
	CHK_EQUAL(b->entries[1].n, 1, errs);
	CHK_CLOSE(b->entries[1].val, 3.0, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);
	return errs;
}

static int testElemOpMul()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;

	/* non-zero a * zero b */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* zero a * non-zero b */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 0, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* non-zero a * non-zero b with zero product */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpAddTo(0, 1, 3.0, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* non-zero a * non-zero b with non-zero product */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, 3.0, &a);
	elemOpAddTo(1, 0, 3.0, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val, 3.0, EPS, errs);
	CHK_EQUAL(b->nOps, 1, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 0, errs);
	CHK_CLOSE(b->entries[0].val, 9.0, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	return errs;
}

static int testSigmaPlus()
{
	int errs = 0;
	ElemOp sp;
	sp = sigmaPlus();
	CHK_EQUAL(sp->nOps, 1, errs);
	CHK_EQUAL(sp->entries[0].m, 1, errs);
	CHK_EQUAL(sp->entries[0].n, 0, errs);
	CHK_CLOSE(sp->entries[0].val, 1.0, EPS, errs);
	elemOpDestroy(&sp);
	return errs;
}

static int testSigmaMinus()
{
	int errs = 0;
	ElemOp sp;
	sp = sigmaMinus();
	CHK_EQUAL(sp->nOps, 1, errs);
	CHK_EQUAL(sp->entries[0].m, 0, errs);
	CHK_EQUAL(sp->entries[0].n, 1, errs);
	CHK_CLOSE(sp->entries[0].val, 1.0, EPS, errs);
	elemOpDestroy(&sp);
	return errs;
}

static int testSigmaZ()
{
	int errs = 0;
	ElemOp sz;
	sz = sigmaZ();
	CHK_EQUAL(sz->nOps, 2, errs);
	CHK_EQUAL(sz->entries[0].m, 1, errs);
	CHK_EQUAL(sz->entries[0].n, 1, errs);
	CHK_CLOSE(sz->entries[0].val, 1.0, EPS, errs);
	CHK_EQUAL(sz->entries[1].m, 0, errs);
	CHK_EQUAL(sz->entries[1].n, 0, errs);
	CHK_CLOSE(sz->entries[1].val, -1.0, EPS, errs);
	elemOpDestroy(&sz);
	return errs;
}

static int testEye()
{
	int errs = 0;
	int i;
	ElemOp e;
	e = eye(5);
	CHK_EQUAL(e->nOps, 5, errs);
	for (i = 0; i < e->nOps; ++i) {
		CHK_EQUAL(e->entries[i].m, i, errs);
		CHK_EQUAL(e->entries[i].n, i, errs);
		CHK_CLOSE(e->entries[i].val, 1.0, EPS, errs);
	}
	elemOpDestroy(&e);
	return errs;
}

static int testNumOp()
{
	int errs = 0;
	int i;
	ElemOp num;
	num = numOp(5);
	CHK_EQUAL(num->nOps, 4, errs);
	for (i = 0; i < num->nOps; ++i) {
		CHK_EQUAL(num->entries[i].m, i + 1, errs);
		CHK_EQUAL(num->entries[i].n, i + 1, errs);
		CHK_CLOSE(num->entries[i].val, i + 1.0, EPS, errs);
	}
	elemOpDestroy(&num);
	return errs;
}

int testAnnihilationOp()
{
	int errs = 0;
	int i;
	ElemOp a;
	a = annihilationOp(5);
	CHK_EQUAL(a->nOps, 4, errs);
	for (i = 0; i < a->nOps; ++i) {
		CHK_EQUAL(a->entries[i].m, i, errs);
		CHK_EQUAL(a->entries[i].n, i + 1, errs);
		CHK_CLOSE(a->entries[i].val, sqrt(i + 1), EPS, errs);
	}
	elemOpDestroy(&a);
	return errs;
}

static int testCreationOp()
{
	int errs = 0;
	int i;
	ElemOp ad;
	ad = creationOp(5);
	CHK_EQUAL(ad->nOps, 4, errs);
	for (i = 0; i < ad->nOps; ++i) {
		CHK_EQUAL(ad->entries[i].m, i + 1, errs);
		CHK_EQUAL(ad->entries[i].n, i, errs);
		CHK_CLOSE(ad->entries[i].val, sqrt(i + 1), EPS, errs);
	}
	elemOpDestroy(&ad);
	return errs;
}

int elemOpTest()
{
	int errs = 0;
	errs += testElemOpCreate();
	errs += testElemOpAddTo();
	errs += testElemOpScale();
	errs += testElemOpPlus();
	errs += testElemOpMul();
	errs += testSigmaPlus();
	errs += testSigmaMinus();
	errs += testSigmaZ();
	errs += testEye();
	errs += testNumOp();
	errs += testAnnihilationOp();
	errs += testCreationOp();
	return errs;
}
