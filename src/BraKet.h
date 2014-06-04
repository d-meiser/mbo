#ifndef BRA_KET_H
#define BRA_KET_H

#include "Amplitude.h"

/*
 * @brief A sparse matrix entry
 *
 */
struct BraKet {
	int m;
	int n;
	Amplitude val;
};

#endif
