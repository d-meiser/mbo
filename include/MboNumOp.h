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
 *
 * MboNumOp objects provide operators in a form suitable for numerical
 * operations such as matrix-vector multiplies.  MboNumOp objects can be
 * obtained from MboTensorOps by means of a compilation step with
 * #mboNumOpCompile.
 *
 * @sa MboNumOp, mboNumOpCompile, MboTensorOp
 *
 * @defgroup mbo_numop MboNumOp
 * @ingroup mbo_core
 * @{
 * */

#ifndef MBO_NUM_OP_H
#define MBO_NUM_OP_H

#include <MboExport.h>
#include <MboTensorOp.h>
#include <MboErrors.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboNumOp_t;
/** @brief Data type for numerical operators. 
 *
 * MboNumOp objects can be created from MboTensorOps using
 * #mboNumOpCompile
 *
 * @sa mboNumOpCompile, MboTensorOp
 * */
typedef struct MboNumOp_t *MboNumOp;

/**
 * @brief Build a numerical operator from a tensor operator.
 *
 * @param op    Operator to compile.
 * @param numOp On exit this pointer contains the compiled operator.
 *              The generated MboNumOp has to be destroyed with
 *              #mboNumOpDestroy.
 *
 * @sa MboNumOp, MboTensorOp, mboNumOpDestroy
 * */
MBO_EXPORT MBO_STATUS mboNumOpCompile(MboTensorOp op, MboNumOp *numOp);

/**
 * @brief Deallocate all resources associated with a MboNumOp.
 *
 * @sa mboNumOpCompile, MboNumOp
 * */
MBO_EXPORT void mboNumOpDestroy(MboNumOp *op);

/**
 * @brief Returns a (non-owning) pointer to the space of an operator.
 *
 * @param op The operator.
 * @return   The product space.
 *
 * @sa MboProdSpace, MboNumOp
 * */
MBO_EXPORT MboProdSpace mboNumOpGetSpace(MboNumOp op);

/**
 * y <- alpha * a * x + beta * y
 * */
MBO_EXPORT MBO_STATUS
mboNumOpMatVec(struct MboAmplitude alpha, MboNumOp a, struct MboAmplitude *x,
	       struct MboAmplitude beta, struct MboAmplitude *y);

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
				   MboGlobInd rmax, MboGlobInd *i);

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
 * @return     Error code.
 *
 * @sa mboNumOpDenseMatrix, mboNumOpRowOffsets.
 * */
MBO_EXPORT MBO_STATUS mboNumOpSparseMatrix(MboNumOp op, MboGlobInd rmin,
					   MboGlobInd rmax, MboGlobInd *i, MboGlobInd *j,
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
/** @} */
#endif
