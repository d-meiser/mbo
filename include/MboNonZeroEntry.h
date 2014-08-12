/**
 * @file MboNonZeroEntry.h
 * @brief Definition of data type for nonzero entries.
 * */
#ifndef MBO_NON_ZERO_ENTRY_H
#define MBO_NON_ZERO_ENTRY_H

#include "MboAmplitude.h"

/*
 * @brief A sparse matrix entry
 *
 */
struct MboNonZeroEntry {
	int m;
	int n;
	struct MboAmplitude val;
};

#endif
