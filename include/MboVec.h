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
/**
 * @file MboVec.h
 * @brief The MboVec API.
 *
 * MboVecs are used to describe vectors in product spaces.
 *
 * @defgroup mbo_vec MboVec
 * @ingroup mbo_utilities
 * @{ 
 * */
#ifndef MBO_VEC_H
#define MBO_VEC_H

#include <MboExport.h>
#include <MboErrors.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude;
struct MboVec_t;
/** @brief Data type for representing vectors */
typedef struct MboVec_t *MboVec;


/** @brief Create MboVec of dimension dim */
MBO_EXPORT MBO_STATUS mboVecCreate(MboGlobInd dim, MboVec *v);

/** @brief Release all resources of a MboVec */
MBO_EXPORT MBO_STATUS mboVecDestroy(MboVec *v);

/** @brief Fill vector with nth unit vector
 * Sets the nth entry in the vector to 1 and all other entries to zero.
 * @param n Dimension in vector to set to 1.
 * @param v The vector to be filled by the unit vector. v must have been
 *          previously created with mboVecCreate.
 * @sa mboVecCreate
 * */
MBO_EXPORT MBO_STATUS mboVecUnitVector(MboGlobInd n, MboVec v);

/** @brief Get dimension of vector */
MBO_EXPORT MboGlobInd mboVecDim(MboVec v);

/** @brief Obtain a read-write view of the vector data
 *
 * * The view has to be released using mboVecReleaseView before the vector can
 * be used again.  Getting a modifiable array may incur additional
 * synchronization overheads.  If the array contents doesn't need to be modified
 * better performance may be obtained with mboVecGetViewR.
 *
 * @see mboVecReleaseView, mboVecGetViewR
 * */
MBO_EXPORT MBO_STATUS mboVecGetViewRW(MboVec v, struct MboAmplitude **array);

/** @brief Obtain a read-only view of vector data
 *
 * The view has to be released using mboVecReleaseView before the vector can be
 * used again.  Any modifications to the array lead to undefined behaviour.  To
 * obtain a modifiable array use mboVecGetViewRW.
 *
 * @see mboVecReleaseView, mboVecGetViewRW
 * */
MBO_EXPORT MBO_STATUS mboVecGetViewR(MboVec v, struct MboAmplitude **array);

/** @brief Release a view of the vector data. */
MBO_EXPORT MBO_STATUS mboVecReleaseView(MboVec v, struct MboAmplitude **array);

/** @brief y <- a * x + y
 * */
MBO_EXPORT MBO_STATUS mboVecAXPY(struct MboAmplitude *a, MboVec x, MboVec y);

/** @brief computes the dot product of two vectors
 * @param x vector on left hand side in dot product ("Bra")
 * @param y vector on right hand size in dot product ("Ket")
 * @param result the dot product
 * */
MBO_EXPORT MBO_STATUS mboVecDot(MboVec x, MboVec y, struct MboAmplitude *result);

/** @brief swap contents of two vectors.
 *
 * The dimensions of the two vectors must match.
 *
 * @param x First vector to swap.
 * @param y Second vector to swap.
 * */
MBO_EXPORT MBO_STATUS mboVecSwap(MboVec x, MboVec y);

/** @brief set vector to a constant */
MBO_EXPORT MBO_STATUS mboVecSet(struct MboAmplitude *a, MboVec x);

/** @brief Fill a vector with random numbers */
MBO_EXPORT MBO_STATUS mboVecSetRandom(MboVec x);

/** @brief Add outer product of vectors
 * @param n    Number of arrays.  n can be obtained from an MboProdSpace
 *             object by means of mboProdSpaceSize.
 * @param dims Length of each array. This array has to be at least of
 *             length n.  dims can be obtained from an MboProdSpace
 *             object by means of mboProdSpaceGetDims.
 * @param vecs Array of arrays of vectors the outer product of which is
 *             taken.  Has to contain at least n
 *             arrays and vecs[i] has to be at least of length
 *             dims[i].
 * @param x    The vector to which the result of the outer product is to
 *             be added.
 * @see mboProdSpaceSize, mboProdSpaceGetDims
 * */
MBO_EXPORT MBO_STATUS
mboVecKron(int n, MboLocInd *dims, struct MboAmplitude **vecs, MboVec x);

/** @brief Apply a function to every entry in a vector
 * @param n    Number of dimensions.
 * @param dims Dimensions of underlying product space.
 * @param f    Function to apply to entries.  As the first two arguments
 *             n and dims are passed to f.  The next integer array is
 *             the set of indices for a given invocation of f.  The last
 *             argument is a function defined context.  The context is
 *             never accessed by mboVecMap.
 * @param ctx  Context for function application.  The context is never
 *             accessed by mboVecMap.
 * @param x    The vector to which the function is applied.
 * */
MBO_EXPORT MBO_STATUS
mboVecMap(int n, MboLocInd *dims,
	  void f(int, MboLocInd *, MboLocInd *, void *, struct MboAmplitude *),
	  void *ctx, MboVec x);

/** @brief Duplicate a vector
 *
 * Generates a copy of x.  The new vector must be destroyed with
 * mboVecDestroy.
 *
 * @param x Vector to copy
 * @param y Upon exit contains the duplicate vector.
 *
 * @sa mboVecDestroy, MboVec
 * */
MBO_EXPORT MBO_STATUS mboVecDuplicate(MboVec x, MboVec *y);

/** @brief Check integrity of MboVec
 * @param v The vector to be checked.
 * @returns The number of errors enountered. 
 *
 * @sa MboVec*/
MBO_EXPORT int mboVecCheck(MboVec v);

#ifdef __cplusplus
}
#endif
/** @} */
#endif
