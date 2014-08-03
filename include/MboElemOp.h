#ifndef MBO_ELEM_OP_H
#define MBO_ELEM_OP_H

#ifdef __cplusplus
extern "C" {
#endif

struct MboElemOp;
typedef struct MboElemOp *MboElemOp;

struct MboAmplitude;

void mboElemOpCreate(MboElemOp *);
void mboElemOpDestroy(MboElemOp *);
void mboElemOpAddTo(int, int, struct MboAmplitude*, MboElemOp *);
void mboElemOpScale(struct MboAmplitude*, MboElemOp);
void mboElemOpPlus(MboElemOp, MboElemOp *);
void mboElemOpMul(MboElemOp, MboElemOp *);
MboElemOp mboElemOpCopy(MboElemOp);
int mboElemOpCheck(MboElemOp);
int mboElemOpTest();

MboElemOp mboSigmaPlus();
MboElemOp mboSigmaMinus();
MboElemOp mboSigmaZ();
MboElemOp mboEye(int);
MboElemOp mboNumOp(int);
MboElemOp mboAnnihilationOp(int);
MboElemOp mboCreationOp(int);

#ifdef __cplusplus
}
#endif

#endif
