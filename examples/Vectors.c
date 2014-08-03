#include <MboVec.h>

int main()
{
	double alpha[2] = {2.0, 3.0};
	long d = 10l;
	int i;
	MboVec x;
	MboVec y;
	mboVecCreate(d, &x);
	mboVecCreate(d, &y);

	for (i = 0; i < 10; ++i) {
		mboVecAXPY(alpha, x, y);
	}

	mboVecDestroy(&x);
	mboVecDestroy(&y);
	return 0;
}

