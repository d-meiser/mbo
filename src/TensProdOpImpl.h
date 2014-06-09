#ifndef TENS_PROD_IMPL_H
#define TENS_PROD_IMPL_H

#include "TensProdOp.h"
#include "Errors.h"

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

struct Embedding *FindEmbedding(int, struct Embedding *);
void MultiplyIthEmbeddings(int i, struct Embedding *,
				  struct Embedding *);
struct Embedding *GatherIthEmbedding(int i, struct Embedding *);

#endif
