#include "ProdSpaceImpl.h"
#include "ProdSpace.h"

#include <stdlib.h>

ProdSpace CreateProdSpace(int d)
{
	ProdSpace sp = (ProdSpace) malloc(sizeof(*sp));
	sp->dim = d;
	sp->next = 0;
	return sp;
}

void DestroyProdSpace(ProdSpace sp)
{
	if (sp->next) {
		DestroyProdSpace(sp->next);
		sp->next = 0;
	}
	free(sp);
}

void MultToProdSpace(ProdSpace a, ProdSpace b)
{
	ProdSpace last = b;
	while (last->next) {
		last = last->next;
	}
	last->next = CopyProdSpace(a);
}

ProdSpace CopyProdSpace(ProdSpace sp)
{
	ProdSpace copy = 0;
	if (sp) {
		copy = (ProdSpace) malloc(sizeof(*copy));
		copy->dim = sp->dim;
		copy->next = CopyProdSpace(sp->next);
	}
	return copy;
}
