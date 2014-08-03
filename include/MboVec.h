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
MBO_STATUS mboVecAXPY(struct MboAmplitude* a, MboVec x, MboVec y);

/** @brief set vector to a constant */
MBO_STATUS mboVecSet(struct MboAmplitude* a, MboVec x);

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
