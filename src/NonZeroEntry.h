#ifndef NON_ZERO_ENTRY_H
#define NON_ZERO_ENTRY_H

#include "Amplitude.h"

/*
 * @brief A sparse matrix entry
 *
 */
struct NonZeroEntry {
	int m;
	int n;
	struct Amplitude val;
};

#endif
