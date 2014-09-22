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
 * @file MboTensorOp.h
 * @brief The MboTensorOp API.
 *
 * MboTensorOps provide functionality for describing tensor product
 * operators including an operator algebra (multiplication, addition,
 * Kronecker products) and application of tensor product operators to
 * vectors.
 *
 * @sa MboVec, mboTensorOpMul, mboTensorOpPlus, mboTensorOpKron,
 * mboTensorOpMatVec
 * */
#ifndef MBO_TENSOR_OP_H
#define MBO_TENSOR_OP_H

#include <MboExport.h>
#include <MboProdSpace.h>
#include <MboElemOp.h>
#include <MboErrors.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboTensorOp_t;
/** Data type for representing tensor operators. */
typedef struct MboTensorOp_t *MboTensorOp;

/** 
 * @brief Create a tensor operator corresponding to the Null operator
 *
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpIdentity
 * */
MBO_EXPORT void mboTensorOpNull(MboProdSpace h, MboTensorOp *top);

/** 
 * @brief Create a tensor operator corresponding to the Identity operator
 *
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpNull
 * */
MBO_EXPORT void mboTensorOpIdentity(MboProdSpace h, MboTensorOp *top);

/** 
 * @brief Destroy a tensor operator object.
 *
 * @sa mboTensorOpCreate, mboTensorOpCreateIdentity
*/
MBO_EXPORT void mboTensorOpDestroy(MboTensorOp *top);

/** 
 * @brief Add an embedding to a tensor operator
 *
 * Takes the elementary operator and embeds it into the product space of
 * the tensor operator at slot i.  Schematically, this can be written as
 * follows:
 *
 *     top += I_0 x I_1 x ... x I_(i-1) x elemop x I_(i+1) x ... x I_N
 *
 * @param elemop Elementary operator to be embedded.
 * @param i      Slot where operator gets embedded into the tensor
 *               operators product space.
 * @param top    Tensor operator to which to add the result.
 * */
MBO_EXPORT void mboTensorOpAddTo(MboElemOp elemop, int i, MboTensorOp top);

/** 
 * @brief Get the space on which the operator is defined 
 *  This method returns a non-owning pointer to the space on which the
 *  tensor operator is defined.  The caller must not modify the product
 *  space or destroy it.
 *
 *  @returns Pointer to product space of operator (non-owning)
 *  */
MBO_EXPORT MboProdSpace mboTensorOpGetSpace(MboTensorOp op);

/** 
 * @brief Adds scaled version of embedded operator to tensor operator.
 * @see mboTensorOpAddTo
 * */
MBO_EXPORT void mboTensorOpAddScaledTo(struct MboAmplitude *alpha,
				    MboElemOp elemop, int i, MboTensorOp top);

/** 
 * @brief Multiply two operators and add result to third
 *      (*c) += a * b
 * */
MBO_EXPORT void mboTensorOpMul(MboTensorOp a, MboTensorOp b, MboTensorOp *c);

/**
 * @brief Add a tensor operator to another operator
 *      (*b) += a
 * */
MBO_EXPORT void mboTensorOpPlus(MboTensorOp a, MboTensorOp *b);

/** 
 * @brief Scale operator
 *
 * Schematically:   *a *= alpha;
 * */
MBO_EXPORT void mboTensorOpScale(struct MboAmplitude *alpha, MboTensorOp *a);

/** 
 * @brief Tensor product of two operators.
 *
 *      *c += ops[0] x ops[1] x ... x ops[n] 
 * */
MBO_EXPORT MBO_STATUS mboTensorOpKron(int n, MboTensorOp *ops, MboTensorOp *c);

/**
 * y <- alpha * a * x + beta * y
 * */
MBO_EXPORT MBO_STATUS
mboTensorOpMatVec(struct MboAmplitude alpha, MboTensorOp a,
		  struct MboAmplitude *x, struct MboAmplitude beta,
		  struct MboAmplitude *y, MboGlobInd rmin,
		  MboGlobInd rmax);

/** 
 * @brief Returns an estimate for the number of floating point operations in a
 *         MatVec.
 *
 * A complex multiply-add is counted as 8 floating point operations.
 *
 * @param a The tensor operator for which the flops are to be computed. */
MBO_EXPORT double mboTensorOpFlops(MboTensorOp a);

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
MBO_EXPORT void mboTensorOpDenseMatrix(MboTensorOp a, struct MboAmplitude *mat);

/**
 * @brief Compute row offsets of sparse matrix.
 *
 * This method computes the row offsets of the sparse matrix representing the
 * tensor operator. This corresponds to the i array in aij or Yale sparse matrix
 * format. We have i[0] = 0, i[1] = nnz[rmin], i[2] = nnz[rmin] + nnz[rmin + 1],
 * etc. In particular, i[rmax - rmin] contains the total number of non-zeros in
 * the row range [rmin, rmax).  Note that this method, like
 * mboTensorOpSparseMatrix, counts duplicate entries separately. The i array is
 *suitable for input to mboTensorOpSparseMatrix.
 *
 * @param op   The tensor.
 * @param rmin Start of row range for which to compute the offsets.
 * @param rmax End of row range for which to compute the offsets.
 * @param i    On exit this array contains the row offsets. The offset of the
 *             first row in the range is set to 0. Thus, in order to obtain the
 *             global offsets one has to add the global offset of the first row
 *             to all offsets in i.
 *
 * @sa mboTensorOpSparseMatrix
 **/
MBO_EXPORT void mboTensorOpRowOffsets(MboTensorOp op, MboGlobInd rmin,
				      MboGlobInd rmax, int *i);

/**
 * @brief Create a sparse matrix representation for a tensor operator.
 *
 * After a call to mboTensorOpRowOffsets and mboTensorOpSparseMatrix the
 * arrays i, j, and a contain the compressed sparse row representation
 * of the matrix for the rows [rmin, rmax).
 *
 * @param op   Create sparse matrix for this operator.
 * @param rmin Start of row range for which to compute the sparse
 *             matrix.
 * @param rmax End of row range.
 * @param i    Offsets array computed by mboTensorOpRowOffsets.
 * @param j    Column indices of nonzero entries.
 * @param a    Non-zero entries.
 *
 * @sa mboTensorOpDenseMatrix, mboTensorOpRowOffsets.
 * */
MBO_EXPORT void mboTensorOpSparseMatrix(MboTensorOp op, MboGlobInd rmin,
					MboGlobInd rmax, int *i, int *j,
					struct MboAmplitude *a);

/**
 * @brief Get diagonal entries from Tensor operator.
 *
 * @param op   The operator.
 * @param rmin First row for which to get the diagonal.
 * @param rmax Last row for which to get the diagonal.
 * @param diag The diagonal. Has to point to memory for at least rmax - rmin
 *             entries.
 *
 * @sa MboTensorOp, mboTensorOpDeleteDiagonal
 * */
MBO_EXPORT void mboTensorOpDiagonal(MboTensorOp op, MboGlobInd rmin,
				    MboGlobInd rmax, struct MboAmplitude *diag);

/**
 * @brief Erase the diagonal from an operator.
 *
 * @param op Operator from which to delete the diagonal.
 *
 * @sa MboTensorOp, mboTensorOpDiagonal
 * */
MBO_EXPORT void mboTensorOpDeleteDiagonal(MboTensorOp op);

/** 
 * @brief Check integrity of tensor operator.
 *
 * @param op The operator to check.
 *
 * @returns the number of errors.
 *
 * @sa MboTensorOp
 * */
MBO_EXPORT int mboTensorOpCheck(MboTensorOp op);

#ifdef __cplusplus
}
#endif
#endif
