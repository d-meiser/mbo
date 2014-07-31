#ifndef TEST_UTILS_H
#define TEST_UTILS_H

double sqrt(double);
double fabs(double);
int printf(const char*, ...);

static void report_error(const char* file, int line)
{
	printf("Error occured at %s(%d).\n", file, line);
}

#define CHK_EQUAL(lhs, rhs, err)                                               \
	do {                                                                   \
		if ((lhs) != (rhs)) {                                          \
			report_error(__FILE__, __LINE__);                      \
			++err;                                                 \
		}                                                              \
	} while (0)

#define CHK_CLOSE(lhs, rhs, eps, err)                                          \
	do {                                                                   \
		if (fabs((lhs) - (rhs)) > eps) {                               \
			report_error(__FILE__, __LINE__);                      \
			++err;                                                 \
		}                                                              \
	} while (0)

#define CHK_NOT_EQUAL(lhs, rhs, err)                                           \
	do {                                                                   \
		if ((lhs) == (rhs)) {                                          \
			report_error(__FILE__, __LINE__);                      \
			++err;                                                 \
		}                                                              \
	} while (0)

#endif
