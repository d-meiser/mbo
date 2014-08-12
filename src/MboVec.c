#include <stdlib.h>
#include <string.h>
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

static void fillWithKron(struct MboAmplitude a, int n, int *dims,
			 struct MboAmplitude **vecs, struct MboAmplitude *arr);

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

MBO_STATUS mboVecDuplicate(MboVec x, MboVec *y)
{
	if (x->mapped) return MBO_VEC_IN_USE;
	MBO_STATUS err = mboVecCreate(mboVecDim(x), y);
	if (err != MBO_SUCCESS) return err;
	memcpy((*y)->array, x->array, mboVecDim(x) * sizeof(*x->array));
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

MBO_STATUS mboVecDot(MboVec x, MboVec y, struct MboAmplitude *result)
{
	struct MboAmplitude *xarr, *yarr;
	int i;
	MBO_STATUS err;

	if (x->dim != y->dim) return MBO_DIMENSIONS_MISMATCH;

	err = mboVecGetViewR(x, &xarr);
	if (err) return err;
	err = mboVecGetViewR(y, &yarr);
	if (err) return err;

	result->re = 0;
	result->im = 0;
	for (i = 0; i < x->dim; ++i) {
		result->re += xarr[i].re * yarr[i].re + xarr[i].im * yarr[i].im;
		result->im += xarr[i].re * yarr[i].im - xarr[i].im * yarr[i].re;
	}

	return MBO_SUCCESS;
}

MBO_STATUS mboVecUnitVector(long n, MboVec x)
{
	long i;
	if (n < 0 || n >= x->dim) return MBO_ILLEGAL_DIMENSION;
	if (x->mapped) return MBO_VEC_IN_USE;
	for (i = 0; i < x->dim; ++i) {
		x->array[i].re = (i == n) ? 1.0 : 0.0;
		x->array[i].im = 0.0;
	}
	return MBO_SUCCESS;
}

MBO_STATUS mboVecSwap(MboVec x, MboVec y)
{
	struct MboAmplitude *tmp;
	if (x->dim != y->dim) return MBO_DIMENSIONS_MISMATCH;
	if (x->mapped || y->mapped) return MBO_VEC_IN_USE;
	tmp = y->array;
	y->array = x->array;
	x->array = tmp;
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

MBO_STATUS mboVecSetRandom(MboVec x)
{
	int i;
	if (x->mapped) return MBO_VEC_IN_USE;
	for (i = 0; i < x->dim; ++i) {
		x->array[i].re = rand() / (double)RAND_MAX;
		x->array[i].im = rand() / (double)RAND_MAX;
	}
	return MBO_SUCCESS;
}

MBO_STATUS mboVecKron(int n, int *dims, struct MboAmplitude **vecs, MboVec x)
{
	int i, dim;
	struct MboAmplitude *tmp, one;

	if (x->mapped) return MBO_VEC_IN_USE;
	dim = 1;
	for (i = 0; i < n; ++i) {
		dim *= dims[i];
	}
	if (dim != x->dim) return MBO_DIMENSIONS_MISMATCH;

	tmp = calloc(dim,  sizeof(*tmp));
	one.re = 1.0;
	one.im = 0.0;
	fillWithKron(one, n, dims, vecs, tmp);
	memcpy(x->array, tmp, dim * sizeof(*tmp));
	free(tmp);
	return MBO_SUCCESS;
}

static void bumpIndices(int n, int *dims, int *indices)
{
	int i = n - 1;
	++indices[i];
	while (indices[i] >= dims[i] && i > 0) {
		indices[i] = 0;
		++indices[i - 1];
		--i;
	}
}

MBO_STATUS mboVecMap(int n, int *dims,
		     void f(int, int *, int *, void *, struct MboAmplitude *),
		     void *ctx, MboVec x)
{
	long i, totalDim;
	int *indices;
	struct MboAmplitude *arr;
	MBO_STATUS err;

	if (n <= 0) return MBO_INVALID_ARGUMENT;
	indices = malloc(n * sizeof(*indices));

	totalDim = 1;
	for (i = 0; i < n; ++i) {
		totalDim *= dims[i];
	}
	err = mboVecGetViewRW(x, &arr);
	if (err) return err;

	for (i = 0; i < n; ++i) {
		indices[i] = 0;
	}

	for (i = 0; i < totalDim; ++i) {
		f(n, dims, indices, ctx, arr + i);
		bumpIndices(n, dims, indices); 
	}

	err = mboVecReleaseView(x, &arr);
	free(indices);

	return err;
}

long mboVecDim(MboVec x)
{
	return x->dim;
}

/* Private helper functions */

static void fillWithKron(struct MboAmplitude a, int n, int *dims,
			 struct MboAmplitude **vecs, struct MboAmplitude *arr)
{
	int blockSize, i;
	struct MboAmplitude b;
	if (n > 1) {
		blockSize = 1;
		for (i = 1; i < n; ++i) {
			blockSize *= dims[i];
		}
		for (i = 0; i < dims[0]; ++i) {
			b.re = a.re * vecs[0][i].re - a.im * vecs[0][i].im;
			b.im = a.re * vecs[0][i].im + a.im * vecs[0][i].re;
			fillWithKron(b, n - 1, dims + 1, vecs + 1, arr);
			arr += blockSize;
		}
	} else {
		for (i = 0; i < dims[0]; ++i) {
			arr[i].re +=
			    a.re * vecs[0][i].re - a.im * vecs[0][i].im;
			arr[i].im +=
			    a.re * vecs[0][i].im + a.im * vecs[0][i].re;
		}
	}
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
	struct MboAmplitude a, *arr;
	MBO_STATUS err;

	a.re = 3.7;
	a.im = 2.0;

	err = mboVecCreate(d, &x);
	err = mboVecCreate(d, &y);
	err = mboVecAXPY(&a, x, y);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	err = mboVecCreate(d, &x);
	err = mboVecCreate(d, &y);
	err = mboVecGetViewR(x, &arr);
	err = mboVecAXPY(&a, x, y);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	mboVecCreate(d + 1, &x);
	err = mboVecCreate(d, &y);
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

static int testMboVecSetRandom()
{
	int errs = 0, d = 2, i;
	MboVec x;
	struct MboAmplitude *arr, zero;
	MBO_STATUS err;

	err = mboVecCreate(d, &x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecGetViewR(x, &arr);
	err = mboVecSetRandom(x);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	err = mboVecReleaseView(x, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecDestroy(&x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);

	err = mboVecCreate(d, &x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	zero.re = 0;
	zero.im = 0;
	err = mboVecSet(&zero, x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecSetRandom(x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecGetViewR(x, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	for (i = 0; i < d; ++i) {
		CHK_TRUE(arr[i].re != 0, errs);
		CHK_TRUE(arr[i].im != 0, errs);
	}
	err = mboVecReleaseView(x, &arr);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	err = mboVecDestroy(&x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);

	return errs;
}

static int testMboVecDuplicate()
{
	int errs = 0, err, i, d = 5;
	MboVec x, y;
	struct MboAmplitude *arr;

	mboVecCreate(d, &x);
	mboVecSetRandom(x);

	mboVecGetViewR(x, &arr);
	err = mboVecDuplicate(x, &y);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	mboVecReleaseView(x, &arr);

	err = mboVecDuplicate(x, &y);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	CHK_EQUAL(x->dim, y->dim, errs);
	for (i = 0; i < x->dim; ++i) {
		CHK_CLOSE(x->array[i].re, y->array[i].re, EPS, errs);
		CHK_CLOSE(x->array[i].im, y->array[i].im, EPS, errs);
	}
	mboVecDestroy(&y);
	mboVecDestroy(&x);

	return errs;
}

static void buildArrays(int n, int *dims, struct MboAmplitude ***arrays)
{
	int i, j;
	*arrays = malloc(n * sizeof(**arrays));
	for (i = 0; i < n; ++i) {
		(*arrays)[i] = malloc(dims[i] * sizeof(***arrays));
		for (j = 0; j < dims[i]; ++j) {
			(*arrays)[i][j].re = 2.0 * i - 3.0 * j;
			(*arrays)[i][j].im = -2.1 * i + 3.1 * j;
		}
	}
}

static void freeArrays(int n, struct MboAmplitude ***arrays)
{
	int i;
	for (i = 0; i < n; ++i) {
		free((*arrays)[i]);
	}
	free(*arrays);
	*arrays = 0;
}

static int testMboVecKron()
{
	int errs = 0, n, dims[10], i, j;
	struct MboAmplitude **arrays, *arr;
	struct MboAmplitude zero, tmp;
	MboVec x;
	MBO_STATUS err;

	zero.re = 0;
	zero.im = 0;

	n = 1;
	dims[0] = 2;
	buildArrays(n, dims, &arrays);
	mboVecCreate(dims[0] + 1, &x);
	err = mboVecKron(n, dims, arrays, x);
	CHK_EQUAL(err, MBO_DIMENSIONS_MISMATCH, errs);
	mboVecDestroy(&x);
	freeArrays(n, &arrays);

	n = 1;
	dims[0] = 2;
	buildArrays(n, dims, &arrays);
	mboVecCreate(dims[0], &x);
	mboVecSet(&zero, x);
	err = mboVecKron(n, dims, arrays, x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	mboVecGetViewR(x, &arr);
	CHK_CLOSE(arr[0].re, arrays[0][0].re, EPS, errs);
	CHK_CLOSE(arr[0].im, arrays[0][0].im, EPS, errs);
	CHK_CLOSE(arr[1].re, arrays[0][1].re, EPS, errs);
	CHK_CLOSE(arr[1].im, arrays[0][1].im, EPS, errs);
	mboVecReleaseView(x, &arr);
	mboVecDestroy(&x);
	freeArrays(n, &arrays);

	n = 2;
	dims[0] = 2;
	dims[1] = 3;
	buildArrays(n, dims, &arrays);
	mboVecCreate(dims[0] * dims[1], &x);
	mboVecSet(&zero, x);
	mboVecKron(n, dims, arrays, x);
	mboVecGetViewR(x, &arr);
	for (i = 0; i < dims[0]; ++i) {
		for (j = 0; j < dims[1]; ++j) {
			tmp.re = arrays[0][i].re * arrays[1][j].re -
				arrays[0][i].im * arrays[1][j].im;
			tmp.im = arrays[0][i].re * arrays[1][j].im +
				arrays[0][i].im * arrays[1][j].re;
			CHK_CLOSE(tmp.re, arr[i * dims[1] + j].re, EPS, errs);
			CHK_CLOSE(tmp.im, arr[i * dims[1] + j].im, EPS, errs);
		}
	}
	mboVecReleaseView(x, &arr);
	mboVecDestroy(&x);
	freeArrays(n, &arrays);

	return errs;
}

static int testMboVecDim()
{
	int errs = 0;
	long d;
	MboVec x;

	mboVecCreate(10l, &x);
	d = mboVecDim(x);
	CHK_EQUAL(d, 10, errs);
	mboVecDestroy(&x);

	return errs;
}

static void f1(int n, int *dims, int *indices, void *ctx,
	       struct MboAmplitude *x)
{
	int i;
	double value = 0.0, b = 1.0;
	for (i = n - 1; i >= 0; --i) {
		value += indices[i] * b;
		b *= 10.0;
	}
	x->re = value;
	x->im = 0;
}

static void f2(int n, int *dims, int *indices, void *ctx,
	       struct MboAmplitude *x)
{
	int *i = (int *)ctx;
	x->re = *i;
	x->im = 0;
	++*i;
}

int testMboVecMap()
{
	int errs = 0, dims[] = {2, 3, 2}, i;
	MboVec x;
	struct MboAmplitude *xarr;

	mboVecCreate(12, &x);
	mboVecMap(3, dims, f1, 0, x);
	mboVecGetViewR(x, &xarr);
	CHK_CLOSE(xarr[0].re, 0, EPS, errs);
	CHK_CLOSE(xarr[0].im, 0, EPS, errs);
	CHK_CLOSE(xarr[1].re, 1, EPS, errs);
	CHK_CLOSE(xarr[2].re, 10, EPS, errs);
	CHK_CLOSE(xarr[3].re, 11, EPS, errs);
	CHK_CLOSE(xarr[4].re, 20, EPS, errs);
	CHK_CLOSE(xarr[5].re, 21, EPS, errs);
	CHK_CLOSE(xarr[6].re, 100, EPS, errs);
	CHK_CLOSE(xarr[7].re, 101, EPS, errs);
	CHK_CLOSE(xarr[8].re, 110, EPS, errs);
	CHK_CLOSE(xarr[9].re, 111, EPS, errs);
	CHK_CLOSE(xarr[10].re, 120, EPS, errs);
	CHK_CLOSE(xarr[11].re, 121, EPS, errs);
	mboVecReleaseView(x, &xarr);
	mboVecDestroy(&x);

	mboVecCreate(12, &x);
	i = 0;
	mboVecMap(3, dims, f2, &i, x);
	mboVecGetViewR(x, &xarr);
	CHK_EQUAL(i, 12, errs);
	for (i = 0; i < 12; ++i) {
		CHK_CLOSE(xarr[i].re, i, EPS, errs);
		CHK_CLOSE(xarr[i].im, 0, EPS, errs);
	}
	mboVecReleaseView(x, &xarr);
	mboVecDestroy(&x);
	return errs;
}

int testBumpIndices()
{
	int errs = 0, indices[3], dims[] = {2, 3, 6};

	indices[0] = 0;
	indices[1] = 0;
	indices[2] = 0;
	bumpIndices(3, dims, indices);
	CHK_EQUAL(indices[0], 0, errs);
	CHK_EQUAL(indices[1], 0, errs);
	CHK_EQUAL(indices[2], 1, errs);

	indices[0] = 0;
	indices[1] = 0;
	indices[2] = 5;
	bumpIndices(3, dims, indices);
	CHK_EQUAL(indices[0], 0, errs);
	CHK_EQUAL(indices[1], 1, errs);
	CHK_EQUAL(indices[2], 0, errs);

	indices[0] = 0;
	indices[1] = 2;
	indices[2] = 5;
	bumpIndices(3, dims, indices);
	CHK_EQUAL(indices[0], 1, errs);
	CHK_EQUAL(indices[1], 0, errs);
	CHK_EQUAL(indices[2], 0, errs);

	indices[0] = 0;
	indices[1] = 2;
	indices[2] = 6;
	bumpIndices(3, dims, indices);
	CHK_EQUAL(indices[0], 1, errs);
	CHK_EQUAL(indices[1], 0, errs);
	CHK_EQUAL(indices[2], 0, errs);

	return errs;
}

int testMboVecDot()
{
	int errs = 0;
	MboVec x, y;
	MBO_STATUS err;
	struct MboAmplitude result, *dummy;

	mboVecCreate(3, &x);
	mboVecCreate(2, &y);
	err = mboVecDot(x, y, &result);
	CHK_EQUAL(err, MBO_DIMENSIONS_MISMATCH, errs);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	mboVecCreate(2, &x);
	mboVecCreate(2, &y);
	mboVecGetViewRW(x, &dummy);
	err = mboVecDot(x, y, &result);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	mboVecCreate(2, &x);
	mboVecGetViewRW(x, &dummy);
	dummy[0].re = 2.0;
	dummy[0].im = 1.0;
	dummy[1].re = 3.0;
	dummy[1].im = -1.0;
	mboVecReleaseView(x, &dummy);
	mboVecCreate(2, &y);
	mboVecGetViewRW(y, &dummy);
	dummy[0].re = 1.0;
	dummy[0].im = 3.0;
	dummy[1].re = 10.0;
	dummy[1].im = 15.0;
	mboVecReleaseView(y, &dummy);
	err = mboVecDot(x, y, &result);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	CHK_CLOSE(result.re, 2.0 * 1.0 + 1.0 * 3.0 + 3.0 * 10.0 + (-1.0) * 15.0,
		  EPS, errs);
	CHK_CLOSE(result.im, 2.0 * 3.0 - 1.0 * 1.0 + 3.0 * 15.0 - (-1.0) * 10.0,
		  EPS, errs);
	printf("result.re == %lf\n", result.re);
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	return errs;
}

int testMboVecUnitVector()
{
	int errs = 0;
	MBO_STATUS err;
	MboVec x;
	struct MboAmplitude *xarr;

	mboVecCreate(2, &x);
	err = mboVecUnitVector(2, x);
	CHK_EQUAL(err, MBO_ILLEGAL_DIMENSION, errs);
	mboVecDestroy(&x);

	mboVecCreate(2, &x);
	mboVecGetViewRW(x, &xarr);
	err = mboVecUnitVector(1, x);
	CHK_EQUAL(err, MBO_VEC_IN_USE, errs);
	mboVecReleaseView(x, &xarr);
	mboVecDestroy(&x);

	mboVecCreate(2, &x);
	err = mboVecUnitVector(1, x);
	CHK_EQUAL(err, MBO_SUCCESS, errs);
	mboVecGetViewRW(x, &xarr);
	CHK_CLOSE(xarr[0].re, 0, EPS, errs);
	CHK_CLOSE(xarr[0].im, 0, EPS, errs);
	CHK_CLOSE(xarr[1].re, 1, EPS, errs);
	CHK_CLOSE(xarr[1].im, 0, EPS, errs);
	mboVecReleaseView(x, &xarr);
	mboVecDestroy(&x);

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
	errs += testMboVecSetRandom();
	errs += testMboVecDuplicate();
	errs += testMboVecKron();
	errs += testMboVecDim();
	errs += testMboVecMap();
	errs += testBumpIndices();
	errs += testMboVecDot();
	errs += testMboVecUnitVector();
	return errs;
}

