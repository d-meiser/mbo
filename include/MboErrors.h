/**
 * @file MboErrors.h
 * @brief Error codes used by the MBO library.
 * */
#ifndef MBO_ERRORS_H
#define MBO_ERRORS_H

enum MBO_STATUS {
	MBO_SUCCESS = 0,
	MBO_OUT_OF_MEMORY,
	MBO_SPACE_MISMATCH,
	MBO_VEC_IN_USE,
	MBO_DIMENSIONS_MISMATCH,
	MBO_INVALID_ARGUMENT,
};
typedef enum MBO_STATUS MBO_STATUS;

#endif
