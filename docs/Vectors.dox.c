/**
@defgroup Vectors Vectors
@ingroup mbo_examples
@{
\par
Illustrates the use of vectors in MBO.

\code
#include <MboVec.h>
#include <MboAmplitude.h>

int main()
{
	struct MboAmplitude alpha = {2.0, 3.0};
	long d = 10l;
	int i;
	MboVec x;
	MboVec y;
	mboVecCreate(d, &x);
	mboVecCreate(d, &y);

	for (i = 0; i < 10; ++i) {
		mboVecAXPY(&alpha, x, y);
	}

	mboVecDestroy(&x);
	mboVecDestroy(&y);
	return 0;
}

\endcode

<h2>Plain source code</h2>


\code

#include <MboVec.h>
#include <MboAmplitude.h>

int main()
{
	struct MboAmplitude alpha = {2.0, 3.0};
	long d = 10l;
	int i;
	MboVec x;
	MboVec y;
	mboVecCreate(d, &x);
	mboVecCreate(d, &y);

	for (i = 0; i < 10; ++i) {
		mboVecAXPY(&alpha, x, y);
	}

	mboVecDestroy(&x);
	mboVecDestroy(&y);
	return 0;
}


\endcode
@}
*/
