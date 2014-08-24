#ifndef UTILITIES_H
#define UTILITIES_H

#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

MboGlobInd computeBlockSize(int N, MboLocInd *dims);

#ifdef __cplusplus
}
#endif

#endif
