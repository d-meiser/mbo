#include "TensProdOp.h"
#include "TensProdOpImpl.h"
#include "TestUtils.h"
#include "ElemOp.h"
#include "ElemOpImpl.h"

#include <stdlib.h>

#define EPS 1.0e-12

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

int test_AddToTOp()
{
	int errs = 0;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	double matrixElement = -3.4;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, matrixElement, &eop);

	CreateTOp(0, &op);
	AddToTOp(eop, 0, op);
	CHK_EQUAL(op->sum->first->i, 0, errs);
	CHK_EQUAL(op->sum->first->op->op.m, 14, errs);
	CHK_EQUAL(op->sum->first->op->op.n, 15, errs);
	CHK_CLOSE(op->sum->first->op->op.val, matrixElement, EPS, errs);
	CHK_EQUAL(op->sum->first->op->next, 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_AddScaledToTOp()
{
	int errs = 0;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	double matrixElement = -3.4;
	double alpha = 2.1;

	h = CreateProdSpace(20);

	CreateElemOp(&eop);
	AddToElemOp(14, 15, -3.4, &eop);

	CreateTOp(0, &op);
	AddScaledToTOp(alpha, eop, 0, op);
	CHK_EQUAL(op->sum->first->i, 0, errs);
	CHK_EQUAL(op->sum->first->op->op.m, 14, errs);
	CHK_EQUAL(op->sum->first->op->op.n, 15, errs);
	CHK_CLOSE(op->sum->first->op->op.val, alpha * matrixElement, EPS, errs);
	CHK_EQUAL(op->sum->first->op->next, 0, errs);
	CHK_EQUAL(op->sum->next, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);
	return errs;
}

int test_FindEmbedding()
{
	int errs = 0;
	int i;
	int N = 10;
	TOp op;
	ElemOp eop;
	ProdSpace h;
	struct Embedding *emb;

	h = CreateProdSpace(20);
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &h);
	}

	CreateElemOp(&eop);
	AddToElemOp(14, 15, -3.4, &eop);
	AddToElemOp(1, 5, 2.0, &eop);
	AddToElemOp(2, 3, 2.0, &eop);

	CreateTOp(h, &op);
	AddToTOp(eop, 0, op);
	AddToTOp(eop, 3, op);

	emb = FindEmbedding(0, op->sum->first);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->first);
	CHK_EQUAL(emb, 0, errs);
	emb = FindEmbedding(3, op->sum->next->first);
	CHK_NOT_EQUAL(emb, 0, errs);
	emb = FindEmbedding(1, op->sum->next->first);
	CHK_EQUAL(emb, 0, errs);

	DestroyTOp(&op);
	DestroyElemOp(&eop);
	DestroyProdSpace(&h);

	return errs;
}

int test_MulTOp()
{
	int errs = 0;
	int i;
	int N = 10;
	TOp op1, op2;
	ElemOp eop1, eop2;
	ProdSpace h;

	h = CreateProdSpace(20);
	for (i = 0; i < N; ++i) {
		MultToProdSpace(h, &h);
	}

	CreateElemOp(&eop1);
	AddToElemOp(14, 15, -3.4, &eop1);
	AddToElemOp(1, 5, 2.0, &eop1);
	AddToElemOp(2, 3, 2.0, &eop1);

	CreateElemOp(&eop2);
	AddToElemOp(8, 3, -3.4, &eop2);
	AddToElemOp(3, 4, 2.0, &eop2);

	CreateTOp(h, &op1);
	AddToTOp(eop1, 0, op1);
	CreateTOp(h, &op2);
	AddToTOp(eop2, 0, op2);

	MulTOp(op1, &op2);

	CHK_EQUAL(op1->sum->first->i, 0, errs);
	CHK_EQUAL(op1->sum->next, 0, errs);
	CHK_EQUAL(op1->sum->first->op->op.m, 2, errs);
	CHK_EQUAL(op1->sum->first->op->op.n, 3, errs);
	CHK_EQUAL(op1->sum->first->op->next->op.m, 1, errs);
	CHK_EQUAL(op1->sum->first->op->next->op.n, 5, errs);
	CHK_EQUAL(op1->sum->first->op->next->next->op.m, 14, errs);
	CHK_EQUAL(op1->sum->first->op->next->next->op.n, 15, errs);
	CHK_EQUAL(op1->sum->first->op->next->next->next, 0, errs);

	CHK_EQUAL(op2->sum->first->i, 0, errs);
	CHK_EQUAL(op2->sum->first->next, 0, errs);
	CHK_EQUAL(op2->sum->first->op->op.m, 2, errs);
	CHK_EQUAL(op2->sum->first->op->op.n, 4, errs);
	CHK_CLOSE(op2->sum->first->op->op.val, 4.0, EPS, errs);

	DestroyTOp(&op1);
	DestroyTOp(&op2);
	DestroyElemOp(&eop1);
	DestroyElemOp(&eop2);
	DestroyProdSpace(&h);
	return errs;
}

int main()
{
	int errs = 0;
	errs += test_CreateTOp();
	errs += test_AddToTOp();
	errs += test_AddScaledToTOp();
	errs += test_FindEmbedding();
//	errs += test_MulTOp();
	return errs;
}
