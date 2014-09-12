#ifndef EMBEDDING_H
#define EMBEDDING_H

#include <MboElemOp.h>
#include <MboAmplitude.h>
#include <Utilities.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Elementary embedding of an operator into a product space.
 * */
struct Embedding
{
	int i;
	MboElemOp op;
};

void destroyEmbedding(struct Embedding *e);
int findEmbedding(int i, int numEmbeddings,
			 struct Embedding *embeddings);
void copyEmbedding(struct Embedding *dest, struct Embedding *src);
/*
 * @brief Search for an embedding.
 *
 * Gathers all embeddings embedded at space i into on slot and returns the index
 * of this slot.  This function reorders entries in embeddings.  Pointers into
 * this array are thus invalidated.
 *
 * @param i The slot to search for.
 * @param numEbeddings Number of embeddings in array.
 * @param embeddings Array of embeddings.
 * @return Index into array the embeddings array where resulting embeddings can
 * be found or a negative number if no suitable embedding has been found.
 * */
int gatherIthEmbedding(int i, int *numEmbeddings,
		       struct Embedding **embeddings);

/**
 * @brief Merge embeddings for each slot
 * Effectively calls gatherIthEmbedding for each index in the *embeddings array.
 * */
void gatherAllEmbeddings(int *numEmbeddings, struct Embedding **embeddings);

/**
 * @brief sort embeddings in ascending order according to i */
void sortEmbeddings(int numEmbeddings, struct Embedding *a);

/**
 * @brief Apply a product of embeddings to a vector */
void applyEmbeddings(int i, int numSpaces, MboLocInd *dims,
			    MboGlobInd blockSize, struct MboAmplitude alpha,
			    int numFactors, struct Embedding *embeddings,
			    struct MboAmplitude *xarr,
			    struct MboAmplitude *yarr);

void embeddingDenseMatrix(int i, int numSpaces, MboLocInd *dims,
			  MboGlobInd blockSize, MboGlobInd dim,
			  struct MboAmplitude alpha, int numFactors,
			  struct Embedding *embeddings,
			  struct MboAmplitude *mat);

void embeddingNonZeros(int i, int numSpaces, MboLocInd *dims,
		       MboGlobInd blockSizeAfter, int numFactors,
		       struct Embedding *embeddings, MboGlobInd rmin,
		       MboGlobInd rmax, MboGlobInd offset, int *nnz);

void embeddingSparseMatrix(int i, int numSpaces, MboLocInd *dims,
			   MboGlobInd blockSizeAfter, struct MboAmplitude alpha,
			   int numFactors, struct Embedding *embeddings, int *I,
			   int *J, struct MboAmplitude *A, int *numInserted);

#ifdef __cplusplus
}
#endif

#endif

