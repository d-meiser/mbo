#include <stdlib.h>
#include <string.h>

#include "MboProdSpace.h"

struct MboProdSpace_t {
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

