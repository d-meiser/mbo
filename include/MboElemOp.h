/**
 * @file MboElemOp.h
 * @brief Definition of elementary operators.
 *
 * Elementary operators are used as building blocks of tensor product
 * operators.  In the context of quantum mechanics they are single
 * particle operators.  Mathematically, they are sparse matrices: A
 * collection of non-zero entries.
 * */
#ifndef MBO_ELEM_OP_H
#define MBO_ELEM_OP_H

#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboElemOp;
typedef struct MboElemOp *MboElemOp;

struct MboAmplitude;

void mboElemOpCreate(MboElemOp *);
void mboElemOpDestroy(MboElemOp *);
void mboElemOpAddTo(MboLocInd, MboLocInd, struct MboAmplitude*, MboElemOp *);
void mboElemOpScale(struct MboAmplitude*, MboElemOp);
void mboElemOpPlus(MboElemOp, MboElemOp *);
void mboElemOpMul(MboElemOp, MboElemOp *);
MboElemOp mboElemOpCopy(MboElemOp);
int mboElemOpNumEntries(MboElemOp);
struct MboNonZeroEntry *mboElemOpGetEntries(MboElemOp);
int mboElemOpCheck(MboElemOp);
int mboElemOpTest();

MboElemOp mboSigmaPlus();
MboElemOp mboSigmaMinus();
MboElemOp mboSigmaZ();
MboElemOp mboEye(MboLocInd);
MboElemOp mboNumOp(MboLocInd);
MboElemOp mboAnnihilationOp(MboLocInd);
MboElemOp mboCreationOp(MboLocInd);

#ifdef __cplusplus
}
#endif

#endif
