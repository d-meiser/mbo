#ifndef MBO_VEC_H
#define MBO_VEC_H

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboVec;
/** @brief Data type for representing vectors */
typedef struct MboVec *MboVec;


/** @brief Create MboVec of dimension dim */
void mboVecCreate(long dim, MboVec *v);

/** @brief Release all resources of a MboVec */
void mboVecDestroy(MboVec *v);

/** @brief Obtain a read-write view of the vector data
 *
 * * The view has to be released using mboVecReleaseView before the vector can
 * be used again.  Getting a modifiable array may incur additional
 * synchronization overheads.  If the array contents doesn't need to be modified
 * better performance may be obtained with mboVecGetViewR.
 *
 * @see mboVecReleaseView, mboVecGetViewR
 * */
int mboVecGetViewRW(MboVec v, struct MboAmplitude **array);

/** @brief Obtain a read-only view of vector data
 *
 * The view has to be released using mboVecReleaseView before the vector can be
 * used again.  Any modifications to the array lead to undefined behaviour.  To
 * obtain a modifiable array use mboVecGetViewRW.
 *
 * @see mboVecReleaseView, mboVecGetViewRW
 * */
int mboVecGetViewR(MboVec v, struct MboAmplitude **array);

/** @brief Release a view of the vector data. */
int mboVecReleaseView(MboVec v, struct MboAmplitude **array);

/** @brief Check integrity of MboVec */
int mboVecCheck();

/** @brief Run MboVec tests */
int mboVecTest();

#ifdef __cplusplus
}
#endif
#endif
