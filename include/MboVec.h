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
void mboVecCreate(int dim, MboVec *v);

/** @brief Release all resources of a MboVec */
void mboVecDestroy(MboVec *v);

/** @brief Run MboVec tests */
int mboVecTest();

#ifdef __cplusplus
}
#endif
#endif
