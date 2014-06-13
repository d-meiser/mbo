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

struct Embedding *FindEmbedding(int i, struct Embedding *list);
void MultiplyEmbeddings(int, struct Embedding *a, struct Embedding *b);

/**
 * @brief A simple product of embeddings
 *
 * A SimpleTOp represents a single term in a sum making up a TOp.
 * */
struct SimpleTOp
{
	struct Embedding *embedding;
	struct SimpleTOp *next;
};

struct Embedding *GatherIthEmbedding(int, struct SimpleTOp *);
void DestroySimpleTOp(struct SimpleTOp *term);

/**
   Data structure for tensor product operators.
   */
struct TOp
{
	/* Product space on which the operator is defined. */
	ProdSpace space;
	/* Pointer to first SimpleTOp in list making up the sum over SimpleTOp
	   operators. */
	struct SimpleTOp *sum;
};


#endif
