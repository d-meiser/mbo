#include "ElemOp.h"
#include "BraKet.h"

int main()
{
	struct BraKet bk;
	bk.m = 1;
	bk.n = 2;
	bk.val = 2.0;
	struct ElemOp *op;
	CreateElemOp(bk, &op);
	DestroyElemOp(op);

	return 0;
}
