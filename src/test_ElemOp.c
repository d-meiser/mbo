#include "ElemOpImpl.h"
#include "BraKet.h"
#include "TestUtils.h"

int testCreateElemOp()
{
	int errs = 0;
	ElemOp a;
	CreateElemOp(&a);

        CHK_EQUAL(0, a, errs);

	DestroyElemOp(&a);
	return errs;
}

int testAddToElemOp()
{
	int errs = 0;
	ElemOp a;
	CreateElemOp(&a);

	AddToElemOp(1, 2, 3.0, &a);

        CHK_EQUAL(a->op.m, 1, errs);
        CHK_EQUAL(a->op.n, 2, errs);
        CHK_EQUAL(a->op.val, 3.0, errs);
        CHK_EQUAL(a->next, 0, errs);

	DestroyElemOp(&a);
	return errs;
}

int main()
{
	int errs = 0;
	errs += testCreateElemOp();
	errs += testAddToElemOp();
	return errs;
}
