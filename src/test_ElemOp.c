#include "ElemOpImpl.h"
#include "BraKet.h"
#include "TestUtils.h"

#define EPS 1.0e-12

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

int testScaleElemOp()
{
	int errs = 0;
	ElemOp a;
	CreateElemOp(&a);

	AddToElemOp(1, 2, 3.0, &a);
        ScaleElemOp(3.4, a);
        CHK_CLOSE(a->op.val, 3.0 * 3.4, EPS, errs);

	DestroyElemOp(&a);
	return errs;
}

int testSigmaPlus()
{
	int errs = 0;
        ElemOp sp;
        sp = sigmaPlus();
        CHK_EQUAL(sp->op.m, 1, errs);
        CHK_EQUAL(sp->op.n, 0, errs);
        CHK_CLOSE(sp->op.val, 1.0, EPS, errs);
        CHK_EQUAL(sp->next, 0, errs);
        DestroyElemOp(&sp);
	return errs;
}

int testSigmaMinus()
{
	int errs = 0;
        ElemOp sp;
        sp = sigmaMinus();
        CHK_EQUAL(sp->op.m, 0, errs);
        CHK_EQUAL(sp->op.n, 1, errs);
        CHK_CLOSE(sp->op.val, 1.0, EPS, errs);
        CHK_EQUAL(sp->next, 0, errs);
        DestroyElemOp(&sp);
	return errs;
}

int testSigmaZ()
{
	int errs = 0;
        ElemOp sz;
        sz = sigmaZ();
        CHK_EQUAL(sz->op.m, 0, errs);
        CHK_EQUAL(sz->op.n, 0, errs);
        CHK_CLOSE(sz->op.val, -1.0, EPS, errs);
        CHK_EQUAL(sz->next->op.m, 1, errs);
        CHK_EQUAL(sz->next->op.n, 1, errs);
        CHK_CLOSE(sz->next->op.val, 1.0, EPS, errs);
        CHK_EQUAL(sz->next->next, 0, errs);
        DestroyElemOp(&sz);
	return errs;
}

int testEye()
{
	int errs = 0;
	int i;
	ElemOp e;
	ElemOp n;
	e = eye(5);
	i = 5;
	for (n = e; n != 0; n = n->next) {
		--i;
		CHK_EQUAL(n->op.m, i, errs);
		CHK_EQUAL(n->op.n, i, errs);
		CHK_CLOSE(n->op.val, 1.0, EPS, errs);
	}
	DestroyElemOp(&e);
	return errs;
}

int testNumOp()
{
	int errs = 0;
	int i;
	ElemOp num;
	ElemOp next;
	num = numOp(5);
	i = 5;
	for (next = num; next != 0; next = next->next) {
		--i;
		CHK_EQUAL(next->op.m, i, errs);
		CHK_EQUAL(next->op.n, i, errs);
		CHK_CLOSE(next->op.val, i, EPS, errs);
	}
	DestroyElemOp(&num);
	return errs;
}

int testAnnihilationOp()
{
	int errs = 0;
	int i;
	ElemOp a;
	ElemOp next;
	a = annihilationOp(5);
	i = 5;
	for (next = a; next != 0; next = next->next) {
		--i;
		CHK_EQUAL(next->op.m, i - 1, errs);
		CHK_EQUAL(next->op.n, i, errs);
		CHK_CLOSE(next->op.val, sqrt(i), EPS, errs);
	}
	DestroyElemOp(&a);
	return errs;
}

int testCreationOp()
{
	int errs = 0;
	int i;
	ElemOp ad;
	ElemOp next;
	ad = creationOp(5);
	i = 5;
	for (next = ad; next != 0; next = next->next) {
		--i;
		CHK_EQUAL(next->op.m, i, errs);
		CHK_EQUAL(next->op.n, i - 1, errs);
		CHK_CLOSE(next->op.val, sqrt(i), EPS, errs);
	}
	DestroyElemOp(&ad);
	return errs;
}

int testPlusElemOp()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;

	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	PlusElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b->op.m, 0, errs);
	CHK_EQUAL(b->op.n, 1, errs);
	CHK_CLOSE(b->op.val, 3.0, EPS, errs);
	CHK_EQUAL(b->next, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &b);
	PlusElemOp(a, &b);
	CHK_EQUAL(a, 0, errs);
	CHK_EQUAL(b->op.m, 0, errs);
	CHK_EQUAL(b->op.n, 1, errs);
	CHK_CLOSE(b->op.val, 3.0, EPS, errs);
	CHK_EQUAL(b->next, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	AddToElemOp(0, 1, 3.0, &b);
	PlusElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b->op.m, 0, errs);
	CHK_EQUAL(b->op.n, 1, errs);
	CHK_CLOSE(b->op.val, 3.0, EPS, errs);
	CHK_EQUAL(b->next->op.m, 0, errs);
	CHK_EQUAL(b->next->op.n, 1, errs);
	CHK_CLOSE(b->next->op.val, 3.0, EPS, errs);
	CHK_EQUAL(b->next->next, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	AddToElemOp(5, 11, -2.0, &b);
	PlusElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b->op.m, 0, errs);
	CHK_EQUAL(b->op.n, 1, errs);
	CHK_CLOSE(b->op.val, 3.0, EPS, errs);
	CHK_EQUAL(b->next->op.m, 5, errs);
	CHK_EQUAL(b->next->op.n, 11, errs);
	CHK_CLOSE(b->next->op.val, -2.0, EPS, errs);
	CHK_EQUAL(b->next->next, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);
	return errs;
}

int testMulElemOp()
{
	int errs = 0;
	ElemOp a = 0;
	ElemOp b = 0;

	/* non-zero a * zero b */
	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	MulElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	/* zero a * non-zero b */
	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &b);
	MulElemOp(a, &b);
	CHK_EQUAL(a, 0, errs);
	CHK_EQUAL(b, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	/* non-zero a * non-zero b with zero product */
	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	AddToElemOp(0, 1, 3.0, &b);
	MulElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	/* non-zero a * non-zero b with non-zero product */
	CreateElemOp(&a);
	CreateElemOp(&b);
	AddToElemOp(0, 1, 3.0, &a);
	AddToElemOp(1, 0, 3.0, &b);
	MulElemOp(a, &b);
	CHK_EQUAL(a->op.m, 0, errs);
	CHK_EQUAL(a->op.n, 1, errs);
	CHK_CLOSE(a->op.val, 3.0, EPS, errs);
	CHK_EQUAL(a->next, 0, errs);
	CHK_EQUAL(b->op.m, 0, errs);
	CHK_EQUAL(b->op.n, 0, errs);
	CHK_CLOSE(b->op.val, 9.0, EPS, errs);
	CHK_EQUAL(b->next, 0, errs);
	DestroyElemOp(&a);
	DestroyElemOp(&b);

	return errs;
}

int main()
{
	int errs = 0;
	errs += testCreateElemOp();
	errs += testAddToElemOp();
	errs += testScaleElemOp();
	errs += testSigmaPlus();
	errs += testSigmaMinus();
	errs += testSigmaZ();
	errs += testEye();
	errs += testNumOp();
	errs += testAnnihilationOp();
	errs += testCreationOp();
	errs += testPlusElemOp();
	errs += testMulElemOp();
	return errs;
}
