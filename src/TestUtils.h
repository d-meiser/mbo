#ifndef TEST_UTILS_H
#define TEST_UTILS_H

double fabs(double);

#define CHK_EQUAL(lhs, rhs, err)                                               \
	do {                                                                   \
		if ((lhs) != (rhs)) ++err;                                     \
	} while (0)

#define CHK_CLOSE(lhs, rhs, eps, err)                                          \
	do {                                                                   \
		if (fabs(lhs - rhs) > eps) ++err;                              \
	} while (0)

#define CHK_NOT_EQUAL(lhs, rhs, err)                                           \
	do {                                                                   \
		if ((lhs) == (rhs)) ++err;                                     \
	} while (0)

#endif
