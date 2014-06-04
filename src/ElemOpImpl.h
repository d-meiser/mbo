#ifndef ELEM_OP_IMPL_H
#define ELEM_OP_IMPL_H

#include "ElemOp.h"
#include "BraKet.h"

struct ElemOp {
	struct BraKet op;
	struct ElemOp *next;
};

#endif
