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

#include <MboSys.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboElemOp;
typedef struct MboElemOp *MboElemOp;

struct MboAmplitude;

MBO_API void mboElemOpCreate(MboElemOp *);
MBO_API void mboElemOpDestroy(MboElemOp *);
MBO_API void mboElemOpAddTo(MboLocInd, MboLocInd, struct MboAmplitude *,
			       MboElemOp *);
MBO_API void mboElemOpScale(struct MboAmplitude*, MboElemOp);
MBO_API void mboElemOpPlus(MboElemOp, MboElemOp *);
MBO_API void mboElemOpMul(MboElemOp, MboElemOp *);
MBO_API MboElemOp mboElemOpCopy(MboElemOp);
MBO_API int mboElemOpNumEntries(MboElemOp);
MBO_API struct MboNonZeroEntry *mboElemOpGetEntries(MboElemOp);
MBO_API int mboElemOpCheck(MboElemOp);
MBO_API int mboElemOpTest();

MBO_API MboElemOp mboSigmaPlus();
MBO_API MboElemOp mboSigmaMinus();
MBO_API MboElemOp mboSigmaZ();
MBO_API MboElemOp mboEye(MboLocInd);
MBO_API MboElemOp mboNumOp(MboLocInd);
MBO_API MboElemOp mboAnnihilationOp(MboLocInd);
MBO_API MboElemOp mboCreationOp(MboLocInd);

#ifdef __cplusplus
}
#endif

#endif
