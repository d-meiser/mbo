#include "TensProdOp.h"
#include "TensProdOpImpl.h"
#include "TestUtils.h"
#include "ElemOp.h"

#include <stdlib.h>

int test_CreateTOp()
{
	int errs = 0;
	TOp op;
	CreateTOp(0, &op);
        CHK_EQUAL(op->space, 0, errs);
        CHK_EQUAL(op->sum, 0, errs);
	DestroyTOp(&op);
	return errs;
}

int test_AddToOp()
{
	int errs = 0;
	TOp op;
        ElemOp eop;

        CreateElemOp(&eop);
        AddToElemOp(14, 15, -3.4, &eop);
	CreateTOp(0, &op);
	DestroyTOp(&op);
        DestroyElemOp(&eop);
	return errs;
}

int main()
{
	int errs = 0;
	errs += test_CreateTOp();
	errs += test_AddToOp();
	return errs;
}
