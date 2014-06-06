#ifndef TENS_PROD_IMPL_H
#define TENS_PROD_IMPL_H

#include "TensProdOp.h"

struct Embedding {
	ElemOp op;
	int i;
	struct Embedding *next;
};

/**
 * Data structure for tensor product operators.
 * */
struct TOp {
	ProdSpace space;
	struct Embedding *sum;
};

#endif
