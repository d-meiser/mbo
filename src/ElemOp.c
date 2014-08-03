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

void elemOpAddTo(int m, int n, struct Amplitude *val, ElemOp *eo)
{
	ElemOp a = *eo;
	a->entries = realloc(a->entries, (a->nOps + 1) * sizeof(*a->entries));
	a->entries[a->nOps].m = m;
	a->entries[a->nOps].n = n;
	a->entries[a->nOps].val.re = val->re;
	a->entries[a->nOps].val.im = val->im;
	++a->nOps;
}

void elemOpScale(struct Amplitude *alpha, ElemOp eo)
{
	int i;
	struct Amplitude tmp;
	for (i = 0; i < eo->nOps; ++i) {
		tmp.re = eo->entries[i].val.re * alpha->re -
			 eo->entries[i].val.im * alpha->im;
		tmp.im = eo->entries[i].val.re * alpha->im +
			 eo->entries[i].val.im * alpha->re;
		eo->entries[i].val = tmp;
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
	struct Amplitude one;

	one.re = 1.0;
	one.im = 0.0;
	elemOpCreate(&sp);
	elemOpAddTo(1, 0, &one, &sp);
	return sp;
}

ElemOp sigmaMinus()
{
	ElemOp sm;
	struct Amplitude one;

	one.re = 1.0;
	one.im = 0.0;
	elemOpCreate(&sm);
	elemOpAddTo(0, 1, &one, &sm);
	return sm;
}

ElemOp sigmaZ()
{
	ElemOp sz;
	struct Amplitude tmp;

	elemOpCreate(&sz);
	tmp.re = 1.0;
	tmp.im = 0.0;
	elemOpAddTo(1, 1, &tmp, &sz);
	tmp.re = -1.0;
	elemOpAddTo(0, 0, &tmp, &sz);
	return sz;
}

ElemOp eye(int d)
{
	ElemOp e;
	int i;
	struct Amplitude tmp;

	elemOpCreate(&e);
	tmp.re = 1.0;
	tmp.im = 0;
	for (i = 0; i < d; ++i) {
		elemOpAddTo(i, i, &tmp, &e);
	}
	return e;
}

ElemOp numOp(int d)
{
	ElemOp n;
	int i;
	struct Amplitude tmp;

	elemOpCreate(&n);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = i;
		elemOpAddTo(i, i, &tmp, &n);
	}
	return n;
}

ElemOp annihilationOp(int d)
{
	ElemOp a;
	int i;
	struct Amplitude tmp;

	elemOpCreate(&a);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = sqrt(i);
		elemOpAddTo(i - 1, i, &tmp, &a);
	}
	return a;
}

ElemOp creationOp(int d)
{
	ElemOp ad;
	int i;
	struct Amplitude tmp;

	elemOpCreate(&ad);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = sqrt(i);
		elemOpAddTo(i, i - 1, &tmp, &ad);
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
				prod->entries[numOps].val.re =
				    a->entries[i].val.re *
					(*b)->entries[j].val.re -
				    a->entries[i].val.im *
					(*b)->entries[j].val.im;
				prod->entries[numOps].val.im =
				    a->entries[i].val.re *
					(*b)->entries[j].val.im +
				    a->entries[i].val.im *
					(*b)->entries[j].val.re;
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
	struct Amplitude alpha;

	elemOpCreate(&a);

	alpha.re = 3.0;
	alpha.im = -5.0;
	elemOpAddTo(1, 2, &alpha, &a);

	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 1, errs);
	CHK_EQUAL(a->entries[0].n, 2, errs);
	CHK_EQUAL(a->entries[0].val.re, alpha.re, errs);
	CHK_EQUAL(a->entries[0].val.im, alpha.im, errs);

	elemOpDestroy(&a);
	return errs;
}

static int testElemOpScale()
{
	int errs = 0;
	ElemOp a;
	struct Amplitude alpha, beta;

	elemOpCreate(&a);

	alpha.re = 3.0;
	alpha.im = 1.0;
	elemOpAddTo(1, 2, &alpha, &a);
	beta.re = -1.0;
	beta.im = 8.0;
	elemOpScale(&beta, a);
	CHK_CLOSE(a->entries[0].val.re, alpha.re * beta.re - alpha.im * beta.im,
		  EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.re * beta.im + alpha.im * beta.re,
		  EPS, errs);

	elemOpDestroy(&a);
	return errs;
}

static int testElemOpPlus()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;
	struct Amplitude alpha, beta;

	elemOpCreate(&a);
	elemOpCreate(&b);
	alpha.re = 11;
	alpha.im = 22;
	elemOpAddTo(0, 1, &alpha, &a);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, alpha.im, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, alpha.im, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &a);
	elemOpAddTo(0, 1, &alpha, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->entries[1].m, 0, errs);
	CHK_EQUAL(b->entries[1].n, 1, errs);
	CHK_CLOSE(b->entries[1].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[1].val.im, alpha.im, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &a);
	beta.re = -2.0;
	beta.im = 100;
	elemOpAddTo(5, 11, &beta, &b);
	elemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 5, errs);
	CHK_EQUAL(b->entries[0].n, 11, errs);
	CHK_CLOSE(b->entries[0].val.re, beta.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, beta.im, EPS, errs);
	CHK_EQUAL(b->entries[1].m, 0, errs);
	CHK_EQUAL(b->entries[1].n, 1, errs);
	CHK_CLOSE(b->entries[1].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[1].val.im, alpha.im, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);
	return errs;
}

static int testElemOpMul()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;
	struct Amplitude alpha;

	alpha.re = 3.0;
	alpha.im = 2.0;

	/* non-zero a * zero b */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &a);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* zero a * non-zero b */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 0, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* non-zero a * non-zero b with zero product */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &a);
	elemOpAddTo(0, 1, &alpha, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	/* non-zero a * non-zero b with non-zero product */
	elemOpCreate(&a);
	elemOpCreate(&b);
	elemOpAddTo(0, 1, &alpha, &a);
	elemOpAddTo(1, 0, &alpha, &b);
	elemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->nOps, 1, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 0, errs);
	CHK_CLOSE(b->entries[0].val.re,
		  alpha.re * alpha.re - alpha.im * alpha.im, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, 2.0 * alpha.re * alpha.im, EPS, errs);
	elemOpDestroy(&a);
	elemOpDestroy(&b);

	return errs;
}

static int testElemOpCheck()
{
	int errs = 0;
	ElemOp o;

	elemOpCreate(&o);
	elemOpCheck(o);
	elemOpDestroy(&o);

	o = sigmaZ();
	elemOpCheck(o);
	elemOpDestroy(&o);

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
	CHK_CLOSE(sp->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sp->entries[0].val.im, 0.0, EPS, errs);
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
	CHK_CLOSE(sp->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sp->entries[0].val.im, 0.0, EPS, errs);
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
	CHK_CLOSE(sz->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sz->entries[0].val.im, 0.0, EPS, errs);
	CHK_EQUAL(sz->entries[1].m, 0, errs);
	CHK_EQUAL(sz->entries[1].n, 0, errs);
	CHK_CLOSE(sz->entries[1].val.re, -1.0, EPS, errs);
	CHK_CLOSE(sz->entries[1].val.im, 0.0, EPS, errs);
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
		CHK_CLOSE(e->entries[i].val.re, 1.0, EPS, errs);
		CHK_CLOSE(e->entries[i].val.im, 0.0, EPS, errs);
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
		CHK_CLOSE(num->entries[i].val.re, i + 1.0, EPS, errs);
		CHK_CLOSE(num->entries[i].val.im, 0.0, EPS, errs);
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
		CHK_CLOSE(a->entries[i].val.re, sqrt(i + 1), EPS, errs);
		CHK_CLOSE(a->entries[i].val.im, 0.0, EPS, errs);
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
		CHK_CLOSE(ad->entries[i].val.re, sqrt(i + 1), EPS, errs);
		CHK_CLOSE(ad->entries[i].val.im, 0.0, EPS, errs);
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
	errs += testElemOpCheck();
	errs += testSigmaPlus();
	errs += testSigmaMinus();
	errs += testSigmaZ();
	errs += testEye();
	errs += testNumOp();
	errs += testAnnihilationOp();
	errs += testCreationOp();
	return errs;
}
