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
#ifndef SIMPLE_T_OPS_H
#define SIMPLE_T_OPS_H

#include <MboAmplitude.h>
#include <MboNonZeroEntry.h>
#include <MboProdSpace.h>
#include <MboErrors.h>
#include <Embedding.h>
#include <Tile.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief A simple product of embeddings
 *
 * A SimpleTOp represents a single term in a sum making up a TOp.
 * */
struct SimpleTOp
{
	int numFactors;
	struct Embedding *embeddings;
};

void destroySimpleTOp(struct SimpleTOp *term);
void multiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);
void kronSimpleTOps(struct SimpleTOp *a, int numSpacesInA, struct SimpleTOp *b,
		    struct SimpleTOp *c);
void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);
void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h,
		    struct SimpleTOp *op);
int checkSimpleTOp(struct SimpleTOp *sa);
MBO_STATUS applySimpleTOpMask(MboProdSpace h, struct MboAmplitude alpha,
			  struct SimpleTOp *a, struct MboAmplitude *x,
			  struct MboAmplitude *y, const struct Tile *mask);
double flopsSimpleTOp(int numSpaces, MboLocInd *dims, struct SimpleTOp *op);

void simpleTOpDenseMatrix(MboProdSpace h, struct SimpleTOp *simpleOp,
			  struct MboAmplitude *mat);
void simpleTOpGetNonZerosPerRow(MboProdSpace h, struct SimpleTOp *simpleOp,
				MboGlobInd rmin, MboGlobInd rmax, int *nnz);
void simpleTOpSparseMatrix(MboProdSpace h, struct SimpleTOp *simpleOp,
			   MboGlobInd rmin, MboGlobInd rmax, int *i, int *j,
			   struct MboAmplitude *a, int *numInserted);
void simpleTOpDiagonal(MboProdSpace h, struct SimpleTOp *simpleOp,
		       MboGlobInd rmin, MboGlobInd rmax,
		       struct MboAmplitude *diag);
void simpleTOpDeleteDiagonal(struct SimpleTOp *simpleOp);
void simpleTOpNormalize(struct SimpleTOp *simpleOp);

#ifdef __cplusplus
}
#endif

#endif
