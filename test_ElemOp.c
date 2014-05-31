#include "ElemOpImpl.h"
#include "BraKet.h"

int testCreateElemOp()
{
	int errs = 0;
	struct BraKet bk;
	bk.m = 1;
	bk.n = 2;
	bk.val = 2.0;
	ElemOp a;

	CreateElemOp(bk, &a);

	errs += (a->op.m != bk.m);
	errs += (a->op.n != bk.n);
	errs += (a->op.val != bk.val);
	errs += (a->next != 0);

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

	errs += (a->op.m != bka.m);
	errs += (a->op.n != bka.n);
	errs += (a->op.val != bka.val);
	errs += (a->next != 0);
	errs += (b->op.m != bkb.m);
	errs += (b->op.n != bkb.n);
	errs += (b->op.val != bkb.val);
	errs += (b->next == 0);
	errs += (b->next->op.m != bka.m);
	errs += (b->next->op.n != bka.n);
	errs += (b->next->op.val != bka.val);
	errs += (b->next->next != 0);

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
