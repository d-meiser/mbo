#include <Quo.h>

int main()
{
	ProdSpace h;
	ProdSpace hTot;
	ElemOp sm;
	ElemOp sp;
	ElemOp sz;
	const int N = 20;
        TOp Hamiltonian;

	CreateElemOp(&sm);
	AddToElemOp(0, 1, 1.0, &sm);
	CreateElemOp(&sp);
	AddToElemOp(1, 0, 1.0, &sp);
	CreateElemOp(&sz);
	AddToElemOp(1, 1, 1.0, &sz);
	AddToElemOp(0, 0, -1.0, &sz);

	h = CreateProdSpace(2);
	hTot = CreateProdSpace(0);
	for (int i = 0; i < N; ++i) {
		MultToProdSpace(h, &hTot);
	}

        CreateTOp(hTot, &Hamiltonian);

	DestroyElemOp(&sm);
	DestroyElemOp(&sp);
	DestroyElemOp(&sz);
	DestroyProdSpace(&h);
	DestroyProdSpace(&hTot);
        DestroyTOp(&Hamiltonian);

	return 0;
}

