/*
Copyright 2014 Dominic Meiser

This file is part of mbo.

mbo is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

mbo is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with mbo.  If not, see <http://www.gnu.org/licenses/>.
*/
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

struct MboVec_t
{
	MboGlobInd dim;
	enum MBO_VEC_MAPPING_STATUS mapped;
	struct MboAmplitude *array;
};

static void fillWithKron(struct MboAmplitude a, int n, MboLocInd *dims,
			 struct MboAmplitude **vecs, struct MboAmplitude *arr);

MBO_STATUS mboVecCreate(MboGlobInd dim, MboVec *v)
{
	if (dim <= 0) return MBO_ILLEGAL_DIMENSION;
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
	if (*v != 0) {
		free((*v)->array);
		free(*v);
	}
	*v = 0;
	return MBO_SUCCESS;
}

MBO_STATUS mboVecDuplicate(MboVec x, MboVec *y)
{
	MBO_STATUS err;
	if (x->mapped) return MBO_VEC_IN_USE;
	err = mboVecCreate(mboVecDim(x), y);
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
	MboGlobInd i;
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

	err = mboVecReleaseView(x, &xarr);
	if (err) return err;
	err = mboVecReleaseView(y, &yarr);
	if (err) return err;

	return MBO_SUCCESS;
}

MBO_STATUS mboVecUnitVector(MboGlobInd n, MboVec x)
{
	MboGlobInd i;
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

MBO_STATUS mboVecKron(int n, MboLocInd *dims, struct MboAmplitude **vecs,
		      MboVec x)
{
	int i;
  MboGlobInd dim;
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

static void bumpIndices(int n, MboLocInd *dims, MboLocInd *indices)
{
	int i = n - 1;
	++indices[i];
	while (indices[i] >= dims[i] && i > 0) {
		indices[i] = 0;
		++indices[i - 1];
		--i;
	}
}

MBO_STATUS mboVecMap(int n, MboLocInd *dims,
		     void f(int, MboLocInd *, MboLocInd *, void *,
			    struct MboAmplitude *),
		     void *ctx, MboVec x)
{
	MboGlobInd i, totalDim;
	MboLocInd *indices;
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

MboGlobInd mboVecDim(MboVec x)
{
	return x->dim;
}

/* Private helper functions */

static void fillWithKron(struct MboAmplitude a, int n, MboLocInd *dims,
			 struct MboAmplitude **vecs, struct MboAmplitude *arr)
{
	MboGlobInd blockSize;
  int i;
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

