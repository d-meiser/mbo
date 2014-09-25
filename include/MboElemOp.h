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
 *
 * @defgroup mbo_elemop MboElemOp
 * @ingroup mbo_core
 * @{ 
 * */
#ifndef MBO_ELEM_OP_H
#define MBO_ELEM_OP_H

#include <MboExport.h>
#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MboElemOp_t;
/**
 * @brief Type for representing elementary operators.
 *
 * MboElemOps are operators acting on an individual subsystem.  Roughly
 * speaking, MboElemOps correspond to single particle operators in a
 * quantum mechanical context.  Mathematically, MboElemOps are
 * represented by an unordered collection of non-zero entries.
 * MboElemOps can be used to created MboTensorOps (many-particle
 * operators) by embedding them into a tensor product space.
 *
 * @sa MboTensorOp, MboProdSpace, MboNonZeroEntry
 * */
typedef struct MboElemOp_t *MboElemOp;

struct MboAmplitude;

/**
 * @brief Create a MboElemOp object.
 *
 * MboElemOps created with mboElemOpCreate have to be destroyed with
 * mboElemOpDestroy.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT void mboElemOpCreate(MboElemOp *elemOp);

/**
 * @brief Destroy a MboElemOp object.
 *
 * @sa mboElemOpCreate, MboElemOp
 * */
MBO_EXPORT void mboElemOpDestroy(MboElemOp *elemOp);

/**
 * @brief Add a non-zero entry to an MboElemOp.
 *
 * @param r      Row index of non-zero entry.
 * @param c      Column index of non-zero entry.
 * @param a      Matrix element of entry
 * @param elemOp Elementary operator matrix to which to add the entry.
 *
 * @sa mboElemOpScale, MboElemOp
 * */
MBO_EXPORT void mboElemOpAddTo(MboLocInd r, MboLocInd c, struct MboAmplitude *a,
			       MboElemOp *elemOp);
/**
 * @brief Rescale an elementary operator.
 *
 * elemOp <- a * elemOp
 *
 * @param a      Factor by which to scale the operator.
 * @param elemOp The operator to scale.
 *
 * @sa mboElemOpPlus, mboElemOpMul, MboElemOp
 * */
MBO_EXPORT void mboElemOpScale(struct MboAmplitude* a, MboElemOp elemOp);

/**
 * @brief Add elementary operators.
 *
 * elemOpB <- elemOpB + elemOpA
 *
 * @param elemOpA Operator to add.
 * @param elemOpB Operator to which to add.
 *
 * @sa mboElemOpScale, mboElemOpMul, MboElemOp
 * */
MBO_EXPORT void mboElemOpPlus(MboElemOp elemOpA, MboElemOp *elemOpB);

/**
 * @brief Add elementary operators.
 *
 * elemOpB <- elemOpA * elemOpB
 *
 * @param elemOpA Operator to multiply with.
 * @param elemOpB Operator to be multiplied.
 *
 * @sa mboElemOpScale, mboElemOpPlus, MboElemOp
 * */
MBO_EXPORT void mboElemOpMul(MboElemOp elemOpA, MboElemOp *elemOpB);

/**
 * @brief copy an elementary operator
 *
 * The copy returned by this function must be destroyed with
 * mboElemOpDestroy.
 *
 * @param elemOp Operator to copy.
 * @return       A copy of elemOp.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboElemOpCopy(MboElemOp elemOp);

/**
 * @brief Get the number of non-zero entries.
 *
 * Note that the operator may contain duplicate entries.
 *
 * @param elemOp Operator of which the number of entries is to be
 *               determined.
 * @return       Number of non-zero entries in elemOp.
 *
 * @sa mboElemOpGetEntries, MboElemOp
 * */
MBO_EXPORT int mboElemOpNumEntries(MboElemOp elemOp);

/**
 * @brief Get the non-zero entries in the operator.
 *
 * It is illegal to add entries or to delete entries or otherwise modify
 * the pointer returned by mboElemOpGetEntries.  Use
 * mboElemOpDeleteEntry to delete an entry. Modifications to the entries
 * in the array are allowed.
 *
 * @param elemOp Get the non-zero entries of this operator.
 * @return       An array with the non-zero entries.
 *
 * @sa mboElemOpNumEntries, MboElemOp, mboElemOpDeleteEntry
 * */
MBO_EXPORT struct MboNonZeroEntry *mboElemOpGetEntries(MboElemOp elemOp);

/**
 * @brief Delete a non-zero entry from an operator.
 *
 * Use mboElemOpNumEntries and mboElemOpGetEntries to determine e.
 *
 * @param elemOp Operator from which to delete an entry.
 * @param e      The index of the entry to be delete.
 *
 * @sa mboElemOpNumEntries, mboElemOpGetEntries, MboElemOp
 * */
MBO_EXPORT void mboElemOpDeleteEntry(MboElemOp elemOp, int e);

/**
 * @brief Check internal integrity of an elementary operator.
 *
 * @return Number of errors found.
 *
 * @sa MboElemOp
 * */
MBO_EXPORT int mboElemOpCheck(MboElemOp elemOp);

/**
 * @brief Spin raising operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @return raising operotr.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboSigmaPlus();

/**
 * @brief Spin lowering operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @return lowering operator.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboSigmaMinus();

/**
 * @brief Sigma z Pauli spin matrix.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @return Sigma-Z matrix
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboSigmaZ();

/**
 * @brief Identity operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @param dim Dimension of space in which to build the identity
 *            operator.
 * @return    The identity operator.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboEye(MboLocInd dim);

/**
 * @brief Quantum mechanical harmonic oscillator number operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @param dim Dimension of space in which to build the number
 *            operator.
 * @return    The number operator.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboNumOp(MboLocInd dim);

/**
 * @brief Quantum mechanical harmonic oscillator annihilation operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @param dim Dimension of space in which to build the annihilation
 *            operator.
 * @return    The annihilation operator.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboAnnihilationOp(MboLocInd dim);

/**
 * @brief Quantum mechanical harmonic oscillator creation operator.
 *
 * The operator returned by this function must be freed with
 * mboElemOpDestroy.
 *
 * @param dim Dimension of space in which to build the creation
 *            operator.
 * @return    The creation operator.
 *
 * @sa mboElemOpDestroy, MboElemOp
 * */
MBO_EXPORT MboElemOp mboCreationOp(MboLocInd dim);

#ifdef __cplusplus
}
#endif

/** @} */
#endif
