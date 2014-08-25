#ifndef SIMPLE_T_OPS_H
#define SIMPLE_T_OPS_H

#include <MboAmplitude.h>
#include <MboNonZeroEntry.h>
#include <MboProdSpace.h>
#include <MboErrors.h>
#include <MboVec.h>
#include <Embedding.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief A simple product of embeddings
 *
 * A SimpleTOp represents a single term in a sum making up a TOp.
 * */
struct SimpleTOp
{
	int numFactors;
	struct Embedding *embeddings;
};

void destroySimpleTOp(struct SimpleTOp *term);
void multiplySimpleTOps(int, struct SimpleTOp *sa, struct SimpleTOp *sb);
void kronSimpleTOps(struct SimpleTOp *a, int numSpacesInA, struct SimpleTOp *b,
		    struct SimpleTOp *c);
void copySimpleTOp(struct SimpleTOp *dest, struct SimpleTOp *src);
void scaleSimpleTOp(struct MboAmplitude *alpha, MboProdSpace h,
		    struct SimpleTOp *op);
int checkSimpleTOp(struct SimpleTOp *sa);
MBO_STATUS applySimpleTOp(MboProdSpace h, struct MboAmplitude *alpha,
			  struct SimpleTOp *a, MboVec x, MboVec y);
double flopsSimpleTOp(int numSpaces, MboLocInd *dims, struct SimpleTOp *op);

#ifdef __cplusplus
}
#endif

#endif
