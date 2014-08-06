#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <MboElemOp.h>
#include <MboNonZeroEntry.h>

struct MboElemOp
{
	int nOps;
	struct MboNonZeroEntry *entries;
};

void mboElemOpCreate(MboElemOp *eo)
{
	*eo = malloc(sizeof(struct MboElemOp));
	(*eo)->nOps = 0;
	(*eo)->entries = 0;
}

void mboElemOpDestroy(MboElemOp *eo)
{
	free((*eo)->entries);
	free(*eo);
	*eo = 0;
}

void mboElemOpAddTo(int m, int n, struct MboAmplitude *val, MboElemOp *eo)
{
	MboElemOp a = *eo;
	a->entries = realloc(a->entries, (a->nOps + 1) * sizeof(*a->entries));
	a->entries[a->nOps].m = m;
	a->entries[a->nOps].n = n;
	a->entries[a->nOps].val.re = val->re;
	a->entries[a->nOps].val.im = val->im;
	++a->nOps;
}

void mboElemOpScale(struct MboAmplitude *alpha, MboElemOp eo)
{
	int i;
	struct MboAmplitude tmp;
	for (i = 0; i < eo->nOps; ++i) {
		tmp.re = eo->entries[i].val.re * alpha->re -
			 eo->entries[i].val.im * alpha->im;
		tmp.im = eo->entries[i].val.re * alpha->im +
			 eo->entries[i].val.im * alpha->re;
		eo->entries[i].val = tmp;
	}
}

void mboElemOpPlus(MboElemOp a, MboElemOp *b)
{
	(*b)->entries = realloc((*b)->entries, ((*b)->nOps + a->nOps) *
						   sizeof(*(*b)->entries));
	memcpy((*b)->entries + (*b)->nOps, a->entries,
	       a->nOps * sizeof(*a->entries));
	(*b)->nOps += a->nOps;
}

int mboElemOpNumEntries(MboElemOp op)
{
	return op->nOps;
}

struct MboNonZeroEntry *mboElemOpGetEntries(MboElemOp op)
{
	return op->entries;
}

int mboElemOpCheck(MboElemOp a)
{
	int errs = 0;
	int i;

	for (i = 0; i < a->nOps; ++i) {
		if (a->entries[i].m < 0) ++errs;
		if (a->entries[i].n < 0) ++errs;
	}
	return errs;
}

MboElemOp mboSigmaPlus()
{
	MboElemOp sp;
	struct MboAmplitude one;

	one.re = 1.0;
	one.im = 0.0;
	mboElemOpCreate(&sp);
	mboElemOpAddTo(1, 0, &one, &sp);
	return sp;
}

MboElemOp mboSigmaMinus()
{
	MboElemOp sm;
	struct MboAmplitude one;

	one.re = 1.0;
	one.im = 0.0;
	mboElemOpCreate(&sm);
	mboElemOpAddTo(0, 1, &one, &sm);
	return sm;
}

MboElemOp mboSigmaZ()
{
	MboElemOp sz;
	struct MboAmplitude tmp;

	mboElemOpCreate(&sz);
	tmp.re = 1.0;
	tmp.im = 0.0;
	mboElemOpAddTo(1, 1, &tmp, &sz);
	tmp.re = -1.0;
	mboElemOpAddTo(0, 0, &tmp, &sz);
	return sz;
}

MboElemOp mboEye(int d)
{
	MboElemOp e;
	int i;
	struct MboAmplitude tmp;

	mboElemOpCreate(&e);
	tmp.re = 1.0;
	tmp.im = 0;
	for (i = 0; i < d; ++i) {
		mboElemOpAddTo(i, i, &tmp, &e);
	}
	return e;
}

MboElemOp mboNumOp(int d)
{
	MboElemOp n;
	int i;
	struct MboAmplitude tmp;

	mboElemOpCreate(&n);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = i;
		mboElemOpAddTo(i, i, &tmp, &n);
	}
	return n;
}

MboElemOp mboAnnihilationOp(int d)
{
	MboElemOp a;
	int i;
	struct MboAmplitude tmp;

	mboElemOpCreate(&a);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = sqrt(i);
		mboElemOpAddTo(i - 1, i, &tmp, &a);
	}
	return a;
}

MboElemOp mboCreationOp(int d)
{
	MboElemOp ad;
	int i;
	struct MboAmplitude tmp;

	mboElemOpCreate(&ad);
	tmp.im = 0;
	for (i = 1; i < d; ++i) {
		tmp.re = sqrt(i);
		mboElemOpAddTo(i, i - 1, &tmp, &ad);
	}
	return ad;
}

void mboElemOpMul(MboElemOp a, MboElemOp *b)
{
	MboElemOp prod;
	int numOps, i, j;

	mboElemOpCreate(&prod);
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
	mboElemOpDestroy(b);
	*b = prod;
}

MboElemOp mboElemOpCopy(MboElemOp a)
{
	MboElemOp copy;
	mboElemOpCreate(&copy);
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

static int testMboElemOpCreate()
{
	int errs = 0;
	MboElemOp a;
	mboElemOpCreate(&a);

	CHK_EQUAL(0, a->nOps, errs);
	CHK_EQUAL(0, a->entries, errs);

	mboElemOpDestroy(&a);
	return errs;
}

static int testMboElemOpAddTo()
{
	int errs = 0;
	MboElemOp a;
	struct MboAmplitude alpha;

	mboElemOpCreate(&a);

	alpha.re = 3.0;
	alpha.im = -5.0;
	mboElemOpAddTo(1, 2, &alpha, &a);

	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 1, errs);
	CHK_EQUAL(a->entries[0].n, 2, errs);
	CHK_EQUAL(a->entries[0].val.re, alpha.re, errs);
	CHK_EQUAL(a->entries[0].val.im, alpha.im, errs);

	mboElemOpDestroy(&a);
	return errs;
}

static int testMboElemOpScale()
{
	int errs = 0;
	MboElemOp a;
	struct MboAmplitude alpha, beta;

	mboElemOpCreate(&a);

	alpha.re = 3.0;
	alpha.im = 1.0;
	mboElemOpAddTo(1, 2, &alpha, &a);
	beta.re = -1.0;
	beta.im = 8.0;
	mboElemOpScale(&beta, a);
	CHK_CLOSE(a->entries[0].val.re, alpha.re * beta.re - alpha.im * beta.im,
		  EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.re * beta.im + alpha.im * beta.re,
		  EPS, errs);

	mboElemOpDestroy(&a);
	return errs;
}

static int testMboElemOpPlus()
{
	int errs = 0;
	MboElemOp a = 0;
	MboElemOp b = 0;
	struct MboAmplitude alpha, beta;

	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	alpha.re = 11;
	alpha.im = 22;
	mboElemOpAddTo(0, 1, &alpha, &a);
	mboElemOpPlus(a, &b);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, alpha.im, EPS, errs);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &b);
	mboElemOpPlus(a, &b);
	CHK_EQUAL(b->entries[0].m, 0, errs);
	CHK_EQUAL(b->entries[0].n, 1, errs);
	CHK_CLOSE(b->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(b->entries[0].val.im, alpha.im, EPS, errs);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &a);
	mboElemOpAddTo(0, 1, &alpha, &b);
	mboElemOpPlus(a, &b);
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
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &a);
	beta.re = -2.0;
	beta.im = 100;
	mboElemOpAddTo(5, 11, &beta, &b);
	mboElemOpPlus(a, &b);
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
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);
	return errs;
}

static int testMboElemOpMul()
{
	int errs = 0;
	MboElemOp a = 0;
	MboElemOp b = 0;
	struct MboAmplitude alpha;

	alpha.re = 3.0;
	alpha.im = 2.0;

	/* non-zero a * zero b */
	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &a);
	mboElemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	/* zero a * non-zero b */
	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &b);
	mboElemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 0, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	/* non-zero a * non-zero b with zero product */
	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &a);
	mboElemOpAddTo(0, 1, &alpha, &b);
	mboElemOpMul(a, &b);
	CHK_EQUAL(a->nOps, 1, errs);
	CHK_EQUAL(a->entries[0].m, 0, errs);
	CHK_EQUAL(a->entries[0].n, 1, errs);
	CHK_CLOSE(a->entries[0].val.re, alpha.re, EPS, errs);
	CHK_CLOSE(a->entries[0].val.im, alpha.im, EPS, errs);
	CHK_EQUAL(b->nOps, 0, errs);
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	/* non-zero a * non-zero b with non-zero product */
	mboElemOpCreate(&a);
	mboElemOpCreate(&b);
	mboElemOpAddTo(0, 1, &alpha, &a);
	mboElemOpAddTo(1, 0, &alpha, &b);
	mboElemOpMul(a, &b);
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
	mboElemOpDestroy(&a);
	mboElemOpDestroy(&b);

	return errs;
}

static int testMboElemOpCheck()
{
	int errs = 0;
	MboElemOp o;

	mboElemOpCreate(&o);
	mboElemOpCheck(o);
	mboElemOpDestroy(&o);

	o = mboSigmaZ();
	mboElemOpCheck(o);
	mboElemOpDestroy(&o);

	return errs;
}

static int testMboSigmaPlus()
{
	int errs = 0;
	MboElemOp sp;

	sp = mboSigmaPlus();
	CHK_EQUAL(sp->nOps, 1, errs);
	CHK_EQUAL(sp->entries[0].m, 1, errs);
	CHK_EQUAL(sp->entries[0].n, 0, errs);
	CHK_CLOSE(sp->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sp->entries[0].val.im, 0.0, EPS, errs);
	mboElemOpDestroy(&sp);
	return errs;
}

static int testMboSigmaMinus()
{
	int errs = 0;
	MboElemOp sp;

	sp = mboSigmaMinus();
	CHK_EQUAL(sp->nOps, 1, errs);
	CHK_EQUAL(sp->entries[0].m, 0, errs);
	CHK_EQUAL(sp->entries[0].n, 1, errs);
	CHK_CLOSE(sp->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sp->entries[0].val.im, 0.0, EPS, errs);
	mboElemOpDestroy(&sp);
	return errs;
}

static int testMboSigmaZ()
{
	int errs = 0;
	MboElemOp sz;

	sz = mboSigmaZ();
	CHK_EQUAL(sz->nOps, 2, errs);
	CHK_EQUAL(sz->entries[0].m, 1, errs);
	CHK_EQUAL(sz->entries[0].n, 1, errs);
	CHK_CLOSE(sz->entries[0].val.re, 1.0, EPS, errs);
	CHK_CLOSE(sz->entries[0].val.im, 0.0, EPS, errs);
	CHK_EQUAL(sz->entries[1].m, 0, errs);
	CHK_EQUAL(sz->entries[1].n, 0, errs);
	CHK_CLOSE(sz->entries[1].val.re, -1.0, EPS, errs);
	CHK_CLOSE(sz->entries[1].val.im, 0.0, EPS, errs);
	mboElemOpDestroy(&sz);
	return errs;
}

static int testMboEye()
{
	int errs = 0;
	int i;
	MboElemOp e;
	e = mboEye(5);
	CHK_EQUAL(e->nOps, 5, errs);
	for (i = 0; i < e->nOps; ++i) {
		CHK_EQUAL(e->entries[i].m, i, errs);
		CHK_EQUAL(e->entries[i].n, i, errs);
		CHK_CLOSE(e->entries[i].val.re, 1.0, EPS, errs);
		CHK_CLOSE(e->entries[i].val.im, 0.0, EPS, errs);
	}
	mboElemOpDestroy(&e);
	return errs;
}

static int testMboNumOp()
{
	int errs = 0;
	int i;
	MboElemOp num;
	num = mboNumOp(5);
	CHK_EQUAL(num->nOps, 4, errs);
	for (i = 0; i < num->nOps; ++i) {
		CHK_EQUAL(num->entries[i].m, i + 1, errs);
		CHK_EQUAL(num->entries[i].n, i + 1, errs);
		CHK_CLOSE(num->entries[i].val.re, i + 1.0, EPS, errs);
		CHK_CLOSE(num->entries[i].val.im, 0.0, EPS, errs);
	}
	mboElemOpDestroy(&num);
	return errs;
}

int testMboAnnihilationOp()
{
	int errs = 0;
	int i;
	MboElemOp a;
	a = mboAnnihilationOp(5);
	CHK_EQUAL(a->nOps, 4, errs);
	for (i = 0; i < a->nOps; ++i) {
		CHK_EQUAL(a->entries[i].m, i, errs);
		CHK_EQUAL(a->entries[i].n, i + 1, errs);
		CHK_CLOSE(a->entries[i].val.re, sqrt(i + 1), EPS, errs);
		CHK_CLOSE(a->entries[i].val.im, 0.0, EPS, errs);
	}
	mboElemOpDestroy(&a);
	return errs;
}

static int testMboCreationOp()
{
	int errs = 0;
	int i;
	MboElemOp ad;
	ad = mboCreationOp(5);
	CHK_EQUAL(ad->nOps, 4, errs);
	for (i = 0; i < ad->nOps; ++i) {
		CHK_EQUAL(ad->entries[i].m, i + 1, errs);
		CHK_EQUAL(ad->entries[i].n, i, errs);
		CHK_CLOSE(ad->entries[i].val.re, sqrt(i + 1), EPS, errs);
		CHK_CLOSE(ad->entries[i].val.im, 0.0, EPS, errs);
	}
	mboElemOpDestroy(&ad);
	return errs;
}

int mboElemOpTest()
{
	int errs = 0;
	errs += testMboElemOpCreate();
	errs += testMboElemOpAddTo();
	errs += testMboElemOpScale();
	errs += testMboElemOpPlus();
	errs += testMboElemOpMul();
	errs += testMboElemOpCheck();
	errs += testMboSigmaPlus();
	errs += testMboSigmaMinus();
	errs += testMboSigmaZ();
	errs += testMboEye();
	errs += testMboNumOp();
	errs += testMboAnnihilationOp();
	errs += testMboCreationOp();
	return errs;
}
