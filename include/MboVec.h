#ifndef MBO_VEC_H
#define MBO_VEC_H

#include <MboErrors.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboVec;
/** @brief Data type for representing vectors */
typedef struct MboVec *MboVec;


/** @brief Create MboVec of dimension dim */
MBO_STATUS mboVecCreate(long dim, MboVec *v);

/** @brief Release all resources of a MboVec */
MBO_STATUS mboVecDestroy(MboVec *v);

/** @brief Get dimension of vector */
long mboVecDim(MboVec v);

/** @brief Obtain a read-write view of the vector data
 *
 * * The view has to be released using mboVecReleaseView before the vector can
 * be used again.  Getting a modifiable array may incur additional
 * synchronization overheads.  If the array contents doesn't need to be modified
 * better performance may be obtained with mboVecGetViewR.
 *
 * @see mboVecReleaseView, mboVecGetViewR
 * */
MBO_STATUS mboVecGetViewRW(MboVec v, struct MboAmplitude **array);

/** @brief Obtain a read-only view of vector data
 *
 * The view has to be released using mboVecReleaseView before the vector can be
 * used again.  Any modifications to the array lead to undefined behaviour.  To
 * obtain a modifiable array use mboVecGetViewRW.
 *
 * @see mboVecReleaseView, mboVecGetViewRW
 * */
MBO_STATUS mboVecGetViewR(MboVec v, struct MboAmplitude **array);

/** @brief Release a view of the vector data. */
MBO_STATUS mboVecReleaseView(MboVec v, struct MboAmplitude **array);

/** @brief y <- a * x + y
 * */
MBO_STATUS mboVecAXPY(struct MboAmplitude *a, MboVec x, MboVec y);

MBO_STATUS mboVecSwap(MboVec x, MboVec y);

/** @brief set vector to a constant */
MBO_STATUS mboVecSet(struct MboAmplitude *a, MboVec x);

MBO_STATUS mboVecSetRandom(MboVec x);

/** @brief Add outer product of vectors
 * @param n    Number of arrays.  n can be obtained from an MboProdSpace
 *             object by means of mboProdSpaceSize.
 * @param dims Length of each array. This array has to be at least of
 *             length n.  dims can be obtained from an MboProdSpace
 *             object by means of mboProdSpaceGetDims.
 * @param vecs Array of arrays of vectors the outer product of which is
 *             taken.  Has to contain at least n
 *             arrays and vecs[i] has to be at least of length
 *             dims[i].
 * @param x    The vector to which the result of the outer product is to
 *             be added.
 * @see mboProdSpaceSize, mboProdSpaceGetDims
 * */
MBO_STATUS mboVecKron(int n, int *dims, struct MboAmplitude **vecs, MboVec x);

MBO_STATUS mboVecDuplicate(MboVec x, MboVec *y);

/** @brief Check integrity of MboVec
 * Returns the number of errors enountered. */
int mboVecCheck();

/** @brief Run MboVec tests
 * Returns the total number of failures encountered. */
int mboVecTest();

#ifdef __cplusplus
}
#endif
#endif
