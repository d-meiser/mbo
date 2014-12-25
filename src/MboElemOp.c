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
#include <math.h>

#include <MboElemOp.h>
#include <MboNonZeroEntry.h>

struct MboElemOp_t
{
	int nOps;
	struct MboNonZeroEntry *entries;
};

MBO_STATUS mboElemOpCreate(MboElemOp *eo)
{
	*eo = malloc(sizeof(**eo));
	if (!*eo) {
		return MBO_OUT_OF_MEMORY;
	}
	(*eo)->nOps = 0;
	(*eo)->entries = 0;
	return MBO_SUCCESS;
}

void mboElemOpDestroy(MboElemOp *eo)
{
	free((*eo)->entries);
	free(*eo);
	*eo = 0;
}

MBO_STATUS mboElemOpAddTo(MboLocInd m, MboLocInd n, struct MboAmplitude *val,
		    MboElemOp *eo)
{
	MboElemOp a = *eo;
	a->entries = realloc(a->entries, (a->nOps + 1) * sizeof(*a->entries));
	if (!a->entries) {
		return MBO_OUT_OF_MEMORY;
	}
	a->entries[a->nOps].m = m;
	a->entries[a->nOps].n = n;
	a->entries[a->nOps].val.re = val->re;
	a->entries[a->nOps].val.im = val->im;
	++a->nOps;
	return MBO_SUCCESS;
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

MBO_STATUS mboElemOpPlus(MboElemOp a, MboElemOp *b)
{
	(*b)->entries = realloc((*b)->entries, ((*b)->nOps + a->nOps) *
						   sizeof(*(*b)->entries));
	if (!(*b)->entries) {
		return MBO_OUT_OF_MEMORY;
	}
	memcpy((*b)->entries + (*b)->nOps, a->entries,
	       a->nOps * sizeof(*a->entries));
	(*b)->nOps += a->nOps;
	return MBO_SUCCESS;
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

MboElemOp mboEye(MboLocInd d)
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

MboElemOp mboNumOp(MboLocInd d)
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

MboElemOp mboAnnihilationOp(MboLocInd d)
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

MboElemOp mboCreationOp(MboLocInd d)
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

MBO_STATUS mboElemOpMul(MboElemOp a, MboElemOp *b)
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
	if (!prod->entries) {
		return MBO_OUT_OF_MEMORY;
	}
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
	return MBO_SUCCESS;
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

void mboElemOpDeleteEntry(MboElemOp op, int e)
{
	if (e < op->nOps) {
		op->entries[e] = op->entries[op->nOps - 1];
		--op->nOps;
		op->entries =
		    realloc(op->entries, op->nOps * sizeof(*op->entries));
	}
}
