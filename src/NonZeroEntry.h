#ifndef NON_ZERO_ENTRY_H
#define NON_ZERO_ENTRY_H

#include "MboAmplitude.h"

/*
 * @brief A sparse matrix entry
 *
 */
struct NonZeroEntry {
	int m;
	int n;
	struct MboAmplitude val;
};

#endif
