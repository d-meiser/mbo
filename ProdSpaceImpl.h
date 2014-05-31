#ifndef PROD_SPACE_IMPL_H
#define PROD_SPACE_IMPL_H

#include "ProdSpace.h"

struct ProdSpace {
	int dim;
	struct ProdSpace *next;
};

#endif				/* PROD_SPACE_H */
