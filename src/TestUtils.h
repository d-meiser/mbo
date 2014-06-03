#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#define CHK_EQUAL(lhs, rhs, err)                                               \
	do {                                                                   \
		if (!((lhs) == (rhs))) ++err;                                  \
	} while (0)

#define CHK_NOT_EQUAL(lhs, rhs, err)                                           \
	do {                                                                   \
		if (!((lhs) != (rhs))) ++err;                                  \
	} while (0)

#endif
