#include <stdio.h>
#include "MboVec.h"
#include "Mbo.h"

static const int numSpins = 3;

int main()
{
	int i, j;
	MboVec x, y;
	MboElemOp sp, sm;
	MboTensorOp Jx;
	MboProdSpace h1, hTot;
	struct MboAmplitude pointFive, one, zero, dotProduct;

	/* Build Hilbert space */
	h1 = mboProdSpaceCreate(2);
	hTot = mboProdSpaceCreate(0);
	for (i = 0; i < numSpins; ++i) {
		mboProdSpaceMul(h1, &hTot);
	}
	mboProdSpaceDestroy(&h1);
	
	/* Build Jx */
	sp = mboSigmaPlus();
	sm = mboSigmaMinus();
	pointFive.re = 0.5;
	pointFive.im = 0.0;
	mboTensorOpNull(hTot, &Jx);
	for (i = 0; i < numSpins; ++i) {
		mboTensorOpAddScaledTo(&pointFive, sm, i, Jx);
		mboTensorOpAddScaledTo(&pointFive, sp, i, Jx);
	}
	mboElemOpDestroy(&sp);
	mboElemOpDestroy(&sm);

	/* Compute matrix */
	one.re = 1.0;
	one.im = 0.0;
	zero.re = 0.0;
	zero.im = 0.0;
	mboVecCreate(mboProdSpaceDim(hTot), &x);
	mboVecCreate(mboProdSpaceDim(hTot), &y);
	for (i = 0; i < mboProdSpaceDim(hTot); ++i) {
		mboVecUnitVector(i, x);
		for (j = 0; j < mboProdSpaceDim(hTot); ++j) {
			mboVecUnitVector(j, y);
			mboTensorOpMatVec(&one, Jx, y, &zero, y); 
			mboVecDot(x, y, &dotProduct);
			printf("%1.1lf %1.1lf  ", dotProduct.re, dotProduct.im);
		}
		printf("\n");
	}
	mboVecDestroy(&x);
	mboVecDestroy(&y);

	/* Deallocate remaining resources */
	mboTensorOpDestroy(&Jx);
	mboProdSpaceDestroy(&hTot);
}
