#ifndef TENS_PROD_IMPL_H
#define TENS_PROD_IMPL_H

#include "TensProdOp.h"
#include "Errors.h"

struct Embedding
{
	ElemOp op;
	int i;
	struct Embedding *next;
};

struct EmbeddingList
{
	struct Embedding *first;
	struct EmbeddingList *next;
};

/**
 * Data structure for tensor product operators.
 * */
struct TOp
{
	ProdSpace space;
	struct EmbeddingList *sum;
};

struct Embedding *FindEmbedding(int, struct Embedding *);
void MultiplyEmbeddings(int, struct Embedding *, struct Embedding *);
struct Embedding *GatherIthEmbedding(int, struct Embedding *);

#endif
