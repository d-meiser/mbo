#ifndef TENSOR_OP_H
#define TENSOR_OP_H

#include "ProdSpace.h"
#include "ElemOp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TensorOp;
/** Data type for representing tensor operators. */
typedef struct TensorOp *TensorOp;

/** @brief Create a tensor operator corresponding to the Null operator
 * The resources used by the object have to be release with
 * tensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see tensorOpDestroy, tensorOpIdentity
 * */
void tensorOpNull(ProdSpace h, TensorOp *top);

/** @brief Create a tensor operator corresponding to the Identity operator
 * The resources used by the object have to be release with
 * tensorOpDestroy.
 * @param h   Underlying product space for the operator.
 * @param top The operator produced.
 * @see tensorOpDestroy, tensorOpNull
 * */
void tensorOpIdentity(ProdSpace h, TensorOp *top);

/** @brief Destroy a tensor operator object.
 * @see tensorOpCreate, tensorOpCreateIdentity
*/
void tensorOpDestroy(TensorOp *top);

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
void tensorOpAddTo(ElemOp elemop, int i, TensorOp top);

/** @brief Adds scaled version of embedded operator to tensor operator.
 * @see tensorOpAddTo
 * */
void tensorOpAddScaledTo(double alpha, ElemOp elemop, int i, TensorOp top);

/** @brief Multiply two operators and add result to third
 *      (*c) += a * b
 * */
void tensorOpMul(TensorOp a, TensorOp b, TensorOp *c);

/** @brief Add a tensor operator to another operator
 *      (*b) += a
 * */
void tensorOpPlus(TensorOp a, TensorOp *b);
void tensorOpScale(double alpha, TensorOp *);
void tensorOpKron(TensorOp, TensorOp, TensorOp *);
int tensorOpCheck(TensorOp);
int tensorOpTest();

#ifdef __cplusplus
}
#endif
#endif
