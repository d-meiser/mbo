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

struct MboProdSpace_t;
/**
 * @brief Type for describing tensor product spaces
 * */
typedef struct MboProdSpace_t *MboProdSpace;

MBO_EXPORT MboProdSpace mboProdSpaceCreate(MboLocInd);
MBO_EXPORT void mboProdSpaceDestroy(MboProdSpace *);
MBO_EXPORT void mboProdSpaceMul(MboProdSpace, MboProdSpace *);
MBO_EXPORT MboProdSpace mboProdSpaceCopy(MboProdSpace);
MBO_EXPORT MboGlobInd mboProdSpaceDim(MboProdSpace);
MBO_EXPORT int mboProdSpaceSize(MboProdSpace);
MBO_EXPORT void mboProdSpaceGetDims(MboProdSpace, int, MboLocInd *);
MBO_EXPORT int mboProdSpaceEqual(MboProdSpace, MboProdSpace);
MBO_EXPORT int mboProdSpaceCheck(MboProdSpace);

#ifdef __cplusplus
}
#endif
#endif
