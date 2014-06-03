#ifndef TENS_PROD_IMPL_H
#define TENS_PROD_IMPL_H

#include "TensProdOp.h"

struct Embedding {
	ElemOp op;
	int i;
        struct Embedding *next;
};

struct TOp {
	ProdSpace space;
        struct Embedding *sum;
};
#endif
