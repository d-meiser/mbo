#ifndef MBO_TENSOR_OP_H
#define MBO_TENSOR_OP_H

#include <MboProdSpace.h>
#include <MboElemOp.h>
#include <MboErrors.h>
#include <MboVec.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboTensorOp;
/** Data type for representing tensor operators. */
typedef struct MboTensorOp *MboTensorOp;

/** @brief Create a tensor operator corresponding to the Null operator
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpIdentity
 * */
void mboTensorOpNull(MboProdSpace h, MboTensorOp *top);

/** @brief Create a tensor operator corresponding to the Identity operator
 * The resources used by the object have to be release with
 * mboTensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see mboTensorOpDestroy, mboTensorOpNull
 * */
void mboTensorOpIdentity(MboProdSpace h, MboTensorOp *top);

/** @brief Destroy a tensor operator object.
 * @see mboTensorOpCreate, mboTensorOpCreateIdentity
*/
void mboTensorOpDestroy(MboTensorOp *top);

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
void mboTensorOpAddTo(MboElemOp elemop, int i, MboTensorOp top);

/** @brief Adds scaled version of embedded operator to tensor operator.
 * @see mboTensorOpAddTo
 * */
void mboTensorOpAddScaledTo(struct MboAmplitude *alpha, MboElemOp elemop, int i,
			    MboTensorOp top);

/** @brief Multiply two operators and add result to third
 *      (*c) += a * b
 *      @todo : Does this allow for aliasing?
 *      (*b) += a * b?
 *      (*a) += a * a?
 *      It should, if it doesn't already.
 * */
void mboTensorOpMul(MboTensorOp a, MboTensorOp b, MboTensorOp *c);

/** @brief Add a tensor operator to another operator
 *      (*b) += a
 * */
void mboTensorOpPlus(MboTensorOp a, MboTensorOp *b);

/** @brief Scale operator
 * Schematically:   *a *= alpha;
 * */
void mboTensorOpScale(struct MboAmplitude *alpha, MboTensorOp *a);

/** @brief Tensor product of two operators.
 *      *c += a x b
 * */
void mboTensorOpKron(MboTensorOp a, MboTensorOp b, MboTensorOp *c);

/**
 * y <- alpha * a * x + beta * y
 * */
MBO_STATUS mboTensorOpMatVec(struct MboAmplitude *alpha, MboTensorOp a,
			     MboVec x, struct MboAmplitude *beta, MboVec y);

/** @brief Check integrity of tensor operator.
 * Returns the number of errors.*/
int mboTensorOpCheck(MboTensorOp);

/** @brief Run tensor operator test suite.
 * Returns the number of errors. */
int mboTensorOpTest();

#ifdef __cplusplus
}
#endif
#endif
