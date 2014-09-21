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
 * @file MboErrors.h
 * @brief Error codes used by the MBO library.
 * */
#ifndef MBO_ERRORS_H
#define MBO_ERRORS_H

/**
 * @brief Error codes used by the MBO library.
 * */
enum MBO_STATUS {
	MBO_SUCCESS = 0,
	MBO_OUT_OF_MEMORY,
	MBO_SPACE_MISMATCH,
	MBO_VEC_IN_USE,
	MBO_DIMENSIONS_MISMATCH,
	MBO_INVALID_ARGUMENT,
	MBO_ILLEGAL_DIMENSION
};
/**
 * @brief Error codes used by the MBO library.
 * */
typedef enum MBO_STATUS MBO_STATUS;

#endif
