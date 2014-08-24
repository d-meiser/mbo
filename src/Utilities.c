#include <Utilities.h>

MboGlobInd computeBlockSize(int N, MboLocInd *dims)
{
	int i;
	MboGlobInd blockSize = 1;
	for (i = 0; i < N; ++i) {
		blockSize *= (MboGlobInd)dims[i];
	}
	return blockSize;
}


