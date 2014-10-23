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
/**
 * @file MboNumSubMatrix.h
 * @brief Submatrices of numerical operators.
 *
 * @sa MboNumOp, mboNumOpCompile, MboTensorOp
 *
 * @defgroup mbo_numoptile MboNumOp
 * @ingroup mbo_core
 * @{
 * */
#ifndef MBO_NUM_SUB_MATRIX_H
#define MBO_NUM_SUB_MATRIX_H

#include <MboExport.h>
#include <MboIndices.h>
#include <MboAmplitude.h>
#include <MboNumOp.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboNumSubMatrix_t;
/** @brief Data type for sub matrices.
 *
 * MboNumSubMatrix objects can be created from MboNumOps using
 * #mboNumSubMatrixCreate.
 *
 * @sa mboNumSubMatrixCreate, MboNumOp.
 *
 * */
typedef struct MboNumSubMatrix_t *MboNumSubMatrix;

/** @brief Create a sub matrix from an existing #MboNumOp
 *
 * @param op   The operator from which to create the submatrix.  This operator
 *             must remain valid until the returned MboNumSubMatrix is destroyed
 *             or until the operator is replaced with #mboNumSubMatrixSetOp.
 * @param rmin First row of submatrix.
 * @param rmax One past last row of submatrix.
 * @param cmin First column of submatrix.
 * @param cmax One past last column of submatrix.
 * @return     The newly create #MboNumSubMatrix.  This must be destroyed with
 *             mboNumSubMatrixDestroy.
 *
 * @sa mboNumSubMatrixDestroy, mboNumSubMatrixSetOp.
 * */
MBO_EXPORT MboNumSubMatrix
    mboNumSubMatrixCreate(MboNumOp op, MboGlobInd rmin, MboGlobInd rmax,
			  MboGlobInd cmin, MboGlobInd cmax);

MBO_EXPORT void mboNumSubMatrixDestroy(MboNumSubMatrix *m);
MBO_EXPORT void mboNumSubMatrixSetTile(MboNumSubMatrix m, MboGlobInd rmin,
				       MboGlobInd rmax, MboGlobInd cmin,
				       MboGlobInd cmax);
MBO_EXPORT MBO_STATUS
    mboNumSubMatrixMatVec(struct MboAmplitude alpha, MboNumSubMatrix m,
			  struct MboAmplitude *x, struct MboAmplitude beta,
			  struct MboAmplitude *y);

#ifdef __cplusplus
}
#endif
/** @} */
#endif

