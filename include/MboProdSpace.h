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
 * @file MboProdSpace.h
 * @brief The MboProdSpace API.
 *
 * MboProdSpaces are used to describe tensor product spaces.
 * */
#ifndef MBO_PROD_SPACE_H
#define MBO_PROD_SPACE_H

#include <MboExport.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Struct for describing tensor product spaces.
 * */
struct MboProdSpace_t;
/**
 * @brief Type for describing tensor product spaces.
 * */
typedef struct MboProdSpace_t *MboProdSpace;

/** 
 * @brief Create a MboProdSpace.
 *
 The created product space has to be destroyed with
 * mboProdSpaceDestroy.
 *
 * @param dim Dimension of the product space.
 * @return    Produc space 
 *
 * @sa MboProdSpace, mboProdSpaceDestroy
 * */
MBO_EXPORT MboProdSpace mboProdSpaceCreate(MboLocInd dim);

/**
 * @brief Destroy a MboProdSpace.
 *
 * @param h Space to be destroyed.
 *
 * @sa MboProdSpace, mboProdSpaceCreate
 * */
MBO_EXPORT void mboProdSpaceDestroy(MboProdSpace *h);

/**
 * @brief Kronecker product of two spaces.
 *
 * h2 <- h1 * h2
 *
 * @param h1 First operator for Kronecker product.
 * @param h2 On input, the second operator for Kronecker product.  On
 *           exit, h2 contains the product of h1 and h2.
 *
 * @sa MboProdspace
 * */
MBO_EXPORT void mboProdSpaceMul(MboProdSpace h1, MboProdSpace *h2);

/**
 * @brief Create a copy of a product space.
 *
 * The returned MboProdSpace must be destroyed with mboProdSpaceDestroy.
 *
 * @param h The product space to copy.
 * @return  The copy of h.
 *
 * @sa mboProdSpaceDestroy, MboProdSpace
 * */
MBO_EXPORT MboProdSpace mboProdSpaceCopy(MboProdSpace h);

/**
 * @brief Dimension of a tensor product space.
 *
 * @param h The dimension of this space is computed.
 * @return  The global dimension of h.
 *
 * @sa MboProdSpace
 * */
MBO_EXPORT MboGlobInd mboProdSpaceDim(MboProdSpace h);

/**
 * @brief Number of factor spaces.
 *
 * @param h The space.
 * @return  Number of spaces. 
 *
 * @sa MboProdSpace
 * */
MBO_EXPORT int mboProdSpaceSize(MboProdSpace h);

/**
 * @brief Returns the dimensions of the elementary subspaces of h
 *
 * @param h    The space.
 * @param n    Number of dimensions desired.  The first n dimensions are
 *             returned in dims.  If mboProdSpaceSize(h) is smaller than n only
 *             mboProdSpaceSize(h) values are written to dims.  
 * @param dims On exit dims contains the first n dimensions.  dims
 *             should point to space for n or mboProdSpaceSize(h)
 *             integers - whichever is smaller.
 * @sa MboProdSpace
 * */
MBO_EXPORT void mboProdSpaceGetDims(MboProdSpace h, int n, MboLocInd *dims);

/**
 * @brief Check product spaces for equality
 *
 * @param h1 First space.
 * @param h2 Second space.
 * @returns  A non-zero value if h1 == h2.
 *
 * @sa MboProdSpace
 * */
MBO_EXPORT int mboProdSpaceEqual(MboProdSpace h1, MboProdSpace h2);

/**
 * @brief Check internal integrety of a product space.
 *
 * @param h
 * @returns The number of errors found in h.
 *
 * @sa MboProdSpace
 * */
MBO_EXPORT int mboProdSpaceCheck(MboProdSpace h);

#ifdef __cplusplus
}
#endif
#endif
