#include <stdlib.h>
#include <MboVec.h>
#include <MboErrors.h>
#include <MboAmplitude.h>

enum MBO_VEC_MAPPING_STATUS {
	MBO_VEC_MAPPING_STATUS_UNMAPPED = 0,
	MBO_VEC_MAPPING_STATUS_MAPPED_R,
	MBO_VEC_MAPPING_STATUS_MAPPED_RW,
};

struct MboVec {
	long dim;
	enum MBO_VEC_MAPPING_STATUS mapped;
	struct MboAmplitude *array;
};

MBO_STATUS mboVecCreate(long dim, MboVec *v)
{
	*v = malloc(sizeof(**v));
	if (*v == 0) return MBO_OUT_OF_MEMORY;
	(*v)->dim = dim;
	(*v)->mapped = MBO_VEC_MAPPING_STATUS_UNMAPPED;
	(*v)->array = malloc(dim * sizeof(*(*v)->array));
	if ((*v)->array == 0) {
		free(*v);
		*v = 0;
		return MBO_OUT_OF_MEMORY;
	}
	return MBO_SUCCESS;
}

MBO_STATUS mboVecDestroy(MboVec *v)
{
	free((*v)->array);
	free(*v);
	*v = 0;
	return MBO_SUCCESS;
}

int mboVecCheck(MboVec v)
{
	int errs = 0;
	if (v->dim < 0) ++errs;
	if (v->mapped > MBO_VEC_MAPPING_STATUS_MAPPED_RW) ++errs;
	return errs;
}

static MBO_STATUS mboVecGetView(MboVec v, struct MboAmplitude **array,
				enum MBO_VEC_MAPPING_STATUS mapping)
{
	if (v->mapped) {
		*array = 0;
		return MBO_VEC_IN_USE;
	}
	*array = v->array;
	v->mapped = mapping;
	return MBO_SUCCESS;
}

MBO_STATUS mboVecGetViewRW(MboVec v, struct MboAmplitude **array)
{
	return mboVecGetView(v, array, MBO_VEC_MAPPING_STATUS_MAPPED_RW);
}

MBO_STATUS mboVecGetViewR(MboVec v, struct MboAmplitude **array)
{
	return mboVecGetView(v, array, MBO_VEC_MAPPING_STATUS_MAPPED_R);
}

MBO_STATUS mboVecReleaseView(MboVec v, struct MboAmplitude **array)
{
	v->mapped = MBO_VEC_MAPPING_STATUS_UNMAPPED;
	*array = 0;
	return MBO_SUCCESS;
}

MBO_STATUS mboVecAXPY(struct MboAmplitude* a, MboVec x, MboVec y)
{
	long i;
	struct MboAmplitude tmp;
	if (x->dim != y->dim) return MBO_DIMENSIONS_MISMATCH;
	if (x->mapped || y->mapped) return MBO_VEC_IN_USE;
	for (i = 0; i < x->dim; ++i) {
		/* By using a temporary this code works also in the case
		 * where x  == y */
		tmp.re = a->re * x->array[i].re - a->im * x->array[i].im;
		tmp.im = a->re * x->array[i].im + a->im * x->array[i].re;
		y->array[i].re += tmp.re;
		y->array[i].im += tmp.im;
	}
	return MBO_SUCCESS;
}

MBO_STATUS mboVecSet(struct MboAmplitude* a, MboVec x)
{
	int i;
	if (x->mapped) return MBO_VEC_IN_USE;
	for (i = 0; i < x->dim; ++i) {
		x->array[i] = *a;
	}
	return MBO_SUCCESS;
}

#include "TestUtils.h"
#define EPS 1.0e-12

static int testMboVecCreate()
{
	int errs = 0;
	MboVec v;

	mboVecCreate(5, &v);
	CHK_TRUE(v != 0, errs);
	mboVecDestroy(&v);

	return errs;
}

static int testMboVecDestroy()
{
	int errs = 0;
	MboVec v;

	mboVecCreate(5, &v);
	mboVecDestroy(&v);
	CHK_EQUAL(v, 0, errs);

	return errs;
}

static int testMboVecCheck()
{
	int errs = 0;
	MboVec v;

	mboVecCreate(5, &v);
	CHK_EQUAL(mboVecCheck(v), 0, errs);
	mboVecDestroy(&v);

	return errs;
}

static int testMboVecGetViewRW()
{
	int errs = 0, err, i, n;
	MboVec v;
	struct MboAmplitude *arr;

	n = 5;
	mboVecCreate(n, &v);
	err = mboVecGetViewRW(v, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	CHK_TRUE(arr != 0, errs);
	for (i = 0; i < n; ++i) {
		arr[i].re = i;
		arr[i].im = -2.0 * i;
	}
	err = mboVecReleaseView(v, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	CHK_TRUE(arr == 0, errs);
	mboVecDestroy(&v);

	return errs;
}

static int testMboVecGetViewR()
{
	int errs = 0, n = 5, err, i;
	MboVec v;
	struct MboAmplitude *arr;

	mboVecCreate(n, &v);
	err = mboVecGetViewR(v, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	CHK_TRUE(arr != 0, errs);
	err = mboVecGetViewR(v, &arr);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	CHK_TRUE(arr == 0, errs);
	err = mboVecReleaseView(v, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);

	err = mboVecGetViewRW(v, &arr);
	for (i = 0; i < n; ++i) {
		arr[i].re = i;
		arr[i].im = -2.0 * i;
	}
	err = mboVecReleaseView(v, &arr);
	err = mboVecGetViewR(v, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	for (i = 0; i < n; ++i) {
		CHK_CLOSE(arr[i].re, i, EPS, errs);
		CHK_CLOSE(arr[i].im, -2.0 * i, EPS, errs);
	}
	mboVecReleaseView(v, &arr);

	mboVecDestroy(&v);

	return errs;
}

static int testMboVecAXPY()
{
	int errs = 0, d = 2;
	MboVec x, y;
	struct MboAmplitude a;
	MBO_STATUS err;

	a.re = 3.7;
	a.im = 2.0;

	err = mboVecCreate(d, &x);
	err = mboVecCreate(d, &y);
	err = mboVecAXPY(&a, x, y);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	mboVecDestroy(&x);
	mboVecCreate(d + 1, &x);
	err = mboVecAXPY(&a, x, y);
	CHK_EQUAL(err, MBO_DIMENSIONS_MISMATCH, errs);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	return errs;
}

static int testMboVecSet()
{
	int errs = 0, d = 2;
	MboVec x;
	struct MboAmplitude a, *arr;
	MBO_STATUS err;

	a.re = 3.7;
	a.im = 2.0;

	err = mboVecCreate(d, &x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecSet(&a, x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecGetViewR(x, &arr);
	err = mboVecSet(&a, x);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	err = mboVecReleaseView(x, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecDestroy(&x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);

	return errs;
}

int mboVecTest()
{
	int errs = 0;
	errs += testMboVecCreate();
	errs += testMboVecDestroy();
	errs += testMboVecCheck();
	errs += testMboVecGetViewRW();
	errs += testMboVecGetViewR();
	errs += testMboVecAXPY();
	errs += testMboVecSet();
	return errs;
}

