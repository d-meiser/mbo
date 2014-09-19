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
 * @file MboElemOp.h
 * @brief Definition of elementary operators.
 *
 * Elementary operators are used as building blocks of tensor product
 * operators.  In the context of quantum mechanics they are single
 * particle operators.  Mathematically, they are sparse matrices: A
 * collection of non-zero entries.
 * */
#ifndef MBO_ELEM_OP_H
#define MBO_ELEM_OP_H

#include <MboExport.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboElemOp_t;
typedef struct MboElemOp_t *MboElemOp;

struct MboAmplitude;

MBO_EXPORT void mboElemOpCreate(MboElemOp *);
MBO_EXPORT void mboElemOpDestroy(MboElemOp *);
MBO_EXPORT void mboElemOpAddTo(MboLocInd, MboLocInd, struct MboAmplitude *,
			       MboElemOp *);
MBO_EXPORT void mboElemOpScale(struct MboAmplitude*, MboElemOp);
MBO_EXPORT void mboElemOpPlus(MboElemOp, MboElemOp *);
MBO_EXPORT void mboElemOpMul(MboElemOp, MboElemOp *);
MBO_EXPORT MboElemOp mboElemOpCopy(MboElemOp);
MBO_EXPORT int mboElemOpNumEntries(MboElemOp);
MBO_EXPORT struct MboNonZeroEntry *mboElemOpGetEntries(MboElemOp);
MBO_EXPORT void mboElemOpDeleteEntry(MboElemOp op, int e);
MBO_EXPORT int mboElemOpCheck(MboElemOp);

MBO_EXPORT MboElemOp mboSigmaPlus();
MBO_EXPORT MboElemOp mboSigmaMinus();
MBO_EXPORT MboElemOp mboSigmaZ();
MBO_EXPORT MboElemOp mboEye(MboLocInd);
MBO_EXPORT MboElemOp mboNumOp(MboLocInd);
MBO_EXPORT MboElemOp mboAnnihilationOp(MboLocInd);
MBO_EXPORT MboElemOp mboCreationOp(MboLocInd);

#ifdef __cplusplus
}
#endif

#endif
