#include <Quo.h>

int main()
{
	ProdSpace h;
	h = CreateProdSpace(2);
	ProdSpace hTot;
	for (int i = 0; i < N; ++i) {
		MultToProdSpace(h, hTot);
	}
	return 0;
}

