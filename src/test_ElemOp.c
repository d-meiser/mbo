#include "ElemOpImpl.h"
#include "BraKet.h"
#include "TestUtils.h"

int testCreateElemOp()
{
	int errs = 0;
	struct BraKet bk;
	bk.m = 1;
	bk.n = 2;
	bk.val = 2.0;
	ElemOp a;

	CreateElemOp(bk, &a);

        CHK_EQUAL(a->op.m, bk.m, errs);
        CHK_EQUAL(a->op.n, bk.n, errs);
        CHK_EQUAL(a->op.val, bk.val, errs);
        CHK_EQUAL(a->next, 0, errs);

	DestroyElemOp(a);
	return errs;
}

int testAddToElemOp()
{
	int errs = 0;
	struct BraKet bka;
	bka.m = 8;
	bka.n = 9;
	bka.val = 4.0;
	ElemOp a;
	CreateElemOp(bka, &a);
	struct BraKet bkb;
	bkb.m = 1;
	bkb.n = 2;
	bkb.val = 2.0;
	ElemOp b;
	CreateElemOp(bkb, &b);

	AddToElemOp(a, b);

        CHK_EQUAL(a->op.m, bka.m, errs);
        CHK_EQUAL(a->op.n, bka.n, errs);
        CHK_EQUAL(a->op.val, bka.val, errs);
        CHK_EQUAL(a->next, 0, errs);
        CHK_EQUAL(b->op.m, bkb.m, errs);
        CHK_EQUAL(b->op.n, bkb.n, errs);
        CHK_EQUAL(b->op.val, bkb.val, errs);
        CHK_NOT_EQUAL(b->next, 0, errs);
        CHK_EQUAL(b->next->op.m, bka.m, errs);
        CHK_EQUAL(b->next->op.n, bka.n, errs);
        CHK_EQUAL(b->next->op.val, bka.val, errs);
        CHK_EQUAL(b->next->next, 0, errs);

	DestroyElemOp(a);
	DestroyElemOp(b);
	return errs;
}

int main()
{
	int errs = 0;
	errs += testCreateElemOp();
	errs += testAddToElemOp();
	return errs;
}
