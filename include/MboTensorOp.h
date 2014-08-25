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
#include <MboVec.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboTensorOp_t;
/** Data type for representing tensor operators. */
typedef struct MboTensorOp_t *MboTensorOp;

/** @brief Create a tensor operator corresponding to the Null operator
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpIdentity
 * */
MBO_EXPORT void mboTensorOpNull(MboProdSpace h, MboTensorOp *top);

/** @brief Create a tensor operator corresponding to the Identity operator
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpNull
 * */
MBO_EXPORT void mboTensorOpIdentity(MboProdSpace h, MboTensorOp *top);

/** @brief Destroy a tensor operator object.
 * @see mboTensorOpCreate, mboTensorOpCreateIdentity
*/
MBO_EXPORT void mboTensorOpDestroy(MboTensorOp *top);

/** @brief Add an embedding to a tensor operator
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

/** @brief Adds scaled version of embedded operator to tensor operator.
 * @see mboTensorOpAddTo
 * */
MBO_EXPORT void mboTensorOpAddScaledTo(struct MboAmplitude *alpha,
				    MboElemOp elemop, int i, MboTensorOp top);

/** @brief Multiply two operators and add result to third
 *      (*c) += a * b
 *      @todo : Does this allow for aliasing?
 *      (*b) += a * b?
 *      (*a) += a * a?
 *      It should, if it doesn't already.
 * */
MBO_EXPORT void mboTensorOpMul(MboTensorOp a, MboTensorOp b, MboTensorOp *c);

/** @brief Add a tensor operator to another operator
 *      (*b) += a
 * */
MBO_EXPORT void mboTensorOpPlus(MboTensorOp a, MboTensorOp *b);

/** @brief Scale operator
 * Schematically:   *a *= alpha;
 * */
MBO_EXPORT void mboTensorOpScale(struct MboAmplitude *alpha, MboTensorOp *a);

/** @brief Tensor product of two operators.
 *      *c += ops[0] x ops[1] x ... x ops[n] 
 * */
MBO_EXPORT MBO_STATUS mboTensorOpKron(int n, MboTensorOp *ops, MboTensorOp *c);

/**
 * y <- alpha * a * x + beta * y
 * */
MBO_EXPORT MBO_STATUS mboTensorOpMatVec(struct MboAmplitude *alpha, MboTensorOp a,
				     MboVec x, struct MboAmplitude *beta,
				     MboVec y);

/** @brief Returns an estimate for the number of floating point operations in a
 *         MatVec.
 *
 * A complex multiply-add is counted as 8 floating point operations.
 *
 * @param a The tensor operator for which the flops are to be computed. */
MBO_EXPORT double mboTensorOpFlops(MboTensorOp a);

/** @brief Dense matrix representation of operator
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

/** @brief Check integrity of tensor operator.
 * Returns the number of errors.*/
MBO_EXPORT int mboTensorOpCheck(MboTensorOp);

#ifdef __cplusplus
}
#endif
#endif
