/**
@defgroup SimpleExample SimpleExample
@ingroup mbo_examples

@{
\par
Basic example: embed spin lowering operator into a many particle
space.

\code
*/
#include <Mbo.h>

int main()
{
	MboElemOp sm = mboSigmaMinus();
	MboProdSpace singleParticleSpace = mboProdSpaceCreate(2);
	MboProdSpace manyParticleSpace = mboProdSpaceCreate(0);
	int N = 20;
	int i;
	for (i = 0; i < N; ++i) {
		mboProdSpaceMul(singleParticleSpace, &manyParticleSpace);
	}

	MboTensorOp smj;
	mboTensorOpNull(manyParticleSpace, &smj);
	int j = 7;
	mboTensorOpAddTo(sm, j, smj);

	mboTensorOpDestroy(&smj);
	mboProdSpaceDestroy(&singleParticleSpace);
	mboProdSpaceDestroy(&manyParticleSpace);
	mboElemOpDestroy(&sm);
}

/**
\endcode
@}
*/

