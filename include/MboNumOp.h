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
 * @file MboNumOp.h
 * @brief Numerical operators.
 * */

#ifndef MBO_NUM_OP_H
#define MBO_NUM_OP_H

#include <MboExport.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboNumOp_t;
/** Data type for numerical operators. */
typedef struct MboNumOp_t *MboNumOp;

/**
 * @brief Build a numerical operator from a tensor operator.
 * */
MBO_EXPORT MboNumOp mboNumOpCompile(MboTensorOp op);
MBO_EXPORT void mboNumOpDestroy(MboNumOp *op);

MBO_EXPORT MboProdSpace mboNumOpGetSpace(MboNumOp op);

/**
 * y <- alpha * a * x + beta * y
 * */
MBO_EXPORT MBO_STATUS
mboNumOpMatVec(struct MboAmplitude alpha, MboNumOp a, struct MboAmplitude *x,
	       struct MboAmplitude beta, struct MboAmplitude *y,
	       MboGlobInd rmin, MboGlobInd rmax);

/**
 * @brief Returns an estimate for the number of floating point operations in a
 *         MatVec.
 *
 * A complex multiply-add is counted as 8 floating point operations.
 *
 * @param a The numerical operator for which the flops are to be computed. */
MBO_EXPORT double mboNumOpFlops(MboNumOp a);

/**
 * @brief Dense matrix representation of operator.
 *
 * This function fills the array mat with the dense matrix
 * representation of the tensor operator.  The matrix is written in row
 * major order
 *
 * @param a   The tensor operator.
 * @param mat User allocated array for the dense matrix.  Has to be at
 *            least of length dim * dim, where dim is the dimension of
 *            the Hilbert space of the operator.
 * */
MBO_EXPORT void mboNumOpDenseMatrix(MboNumOp a, struct MboAmplitude *mat);

/**
 * @brief Compute row offsets of sparse matrix.
 *
 * This method computes the row offsets of the sparse matrix representing the
 * tensor operator. This corresponds to the i array in aij or Yale sparse matrix
 * format. We have i[0] = 0, i[1] = nnz[rmin], i[2] = nnz[rmin] + nnz[rmin + 1],
 * etc. In particular, i[rmax - rmin] contains the total number of non-zeros in
 * the row range [rmin, rmax).  Note that this method, like
 * mboNumOpSparseMatrix, counts duplicate entries separately. The i array is
 * suitable for input to mboNumOpSparseMatrix.
 *
 * @param op   The tensor.
 * @param rmin Start of row range for which to compute the offsets.
 * @param rmax End of row range for which to compute the offsets.
 * @param i    On exit this array contains the row offsets. The offset of the
 *             first row in the range is set to 0. Thus, in order to obtain the
 *             global offsets one has to add the global offset of the first row
 *             to all offsets in i.
 *
 * @sa mboNumOpSparseMatrix
 **/
MBO_EXPORT void mboNumOpRowOffsets(MboNumOp op, MboGlobInd rmin,
				   MboGlobInd rmax, int *i);

/**
 * @brief Create a sparse matrix representation for a numerical operator.
 *
 * After a call to mboNumOpRowOffsets and mboNumOpSparseMatrix the
 * arrays i, j, and a contain the compressed sparse row representation
 * of the matrix for the rows [rmin, rmax).
 *
 * @param op   Create sparse matrix for this operator.
 * @param rmin Start of row range for which to compute the sparse
 *             matrix.
 * @param rmax End of row range.
 * @param i    Offsets array computed by mboNumOpRowOffsets.
 * @param j    Column indices of nonzero entries.
 * @param a    Non-zero entries.
 *
 * @sa mboNumOpDenseMatrix, mboNumOpRowOffsets.
 * */
MBO_EXPORT void mboNumOpSparseMatrix(MboNumOp op, MboGlobInd rmin,
				     MboGlobInd rmax, int *i, int *j,
				     struct MboAmplitude *a);

/**
 * @brief Get diagonal entries from a numerical operator.
 *
 * @param op   The operator.
 * @param rmin First row for which to get the diagonal.
 * @param rmax Last row for which to get the diagonal.
 * @param diag The diagonal. Has to point to memory for at least rmax - rmin
 *             entries.
 *
 * @sa MboNumOp, mboNumOpDeleteDiagonal
 * */
MBO_EXPORT void mboNumOpDiagonal(MboNumOp op, MboGlobInd rmin, MboGlobInd rmax,
				 struct MboAmplitude *diag);

/**
 * @brief Erase the diagonal from an operator.
 *
 * @param op Operator from which to delete the diagonal.
 *
 * @sa MboNumOp, mboTensorOpDiagonal
 * */
MBO_EXPORT void mboNumOpDeleteDiagonal(MboNumOp op);

#ifdef __cplusplus
}
#endif

#endif
