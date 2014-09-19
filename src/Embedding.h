/*
Copyright 2014 Dominic Meiser

This file is part of mbo.

mbo is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

mbo is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with mbo.  If not, see <http://www.gnu.org/licenses/>.
*/
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
		     struct MboAmplitude *x,
		     struct MboAmplitude *y, MboGlobInd rmin,
		     MboGlobInd rmax);

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
			   int numFactors, struct Embedding *embeddings,
			   MboGlobInd rmin, MboGlobInd rmax, int *I, int *J,
			   struct MboAmplitude *A, MboGlobInd offsetR,
			   MboGlobInd offsetC, int *numInserted);

void embeddingDiagonal(int i, int numSpaces, MboLocInd *dims,
		       MboGlobInd blockSizeAfter, struct MboAmplitude alphak,
		       int numFactors, struct Embedding *embeddings,
		       MboGlobInd rmin, MboGlobInd rmax,
		       struct MboAmplitude *diag, MboGlobInd offsetR,
		       MboGlobInd offsetC);

void embeddingDeleteDiagonal(struct Embedding *embedding);

#ifdef __cplusplus
}
#endif

#endif

