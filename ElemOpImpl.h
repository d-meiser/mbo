#ifndef ELEM_OP_IMPL_H
#define ELEM_OP_IMPL_H

#include "ElemOp.h"

struct ElemOp {
	struct BraKet op;
	struct ElemOp *next;
};

ElemOp ElemOpLast(ElemOp);
ElemOp ElemOpCopy(ElemOp);

#endif
