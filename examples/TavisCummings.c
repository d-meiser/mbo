#include <Quo.h>

int main()
{
	ProdSpace h;
	ProdSpace hTot;
	ElemOp sm;
	ElemOp sp;
	ElemOp sz;
        struct BraKet bk;

	h = CreateProdSpace(2);
	hTot = CreateProdSpace(0);
	for (int i = 0; i < N; ++i) {
		MultToProdSpace(h, hTot);
	}

        CreateElemOp(&sm);
        AddToElemOp(0, 1, 1.0, sm);
        CreateElemOp(&sp);
        AddToElemOp(1, 0, 1.0, sp);
        CreateElemOp(bk, &sz);
        AddToElemOp(1, 1, 1,0, sz);
        AddToElemOp(0, 0, -1.0, sz);

        DestroyElemOp(&sm);
        DestroyElemOp(&sp);
        DestroyElemOp(&sz);
        Destroy(h);
        Destroy(hTot);

	return 0;
}

