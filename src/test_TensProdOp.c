#include "TensProdOp.h"
#include "TensProdOpImpl.h"
#include "TestUtils.h"

#include <stdlib.h>

int test_CreateTOp()
{
	int errs = 0;
	TOp op = CreateTOp(0);
	DestroyTOp(op);
	return errs;
}

int test_AddToOp()
{
	int errs = 0;
	TOp op = CreateTOp(0);
	DestroyTOp(op);
	return errs;
}

int main()
{
	int errs = 0;
	errs += test_CreateTOp();
	errs += test_AddToOp();
	return errs;
}
