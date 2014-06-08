#include "ProdSpaceImpl.h"
#include "ProdSpace.h"

#include <stdlib.h>

ProdSpace CreateProdSpace(int d)
{
	ProdSpace sp = malloc(sizeof(*sp));
	sp->dim = d;
	sp->next = 0;
	return sp;
}

void DestroyProdSpace(ProdSpace *sp)
{
	if ((*sp)->next) {
		DestroyProdSpace(&(*sp)->next);
	}
	free(*sp);
        *sp = 0;
}

void MultToProdSpace(ProdSpace a, ProdSpace *b)
{
	if (a->dim == 0) return;
	if ((*b)->dim == 0) {
		free(*b);
		*b = 0;
	}
	ProdSpace aCopy = CopyProdSpace(a);
        ProdSpace last = aCopy;
        while (last->next) {
		last = last->next;
	}
        last->next = *b;
        *b = aCopy;
}

ProdSpace CopyProdSpace(ProdSpace sp)
{
	ProdSpace copy = 0;
	if (sp) {
		copy = malloc(sizeof(*copy));
		copy->dim = sp->dim;
		copy->next = CopyProdSpace(sp->next);
	}
	return copy;
}

long long DimProdSpace(ProdSpace sp)
{
	long long dim = 1;
	while (sp) {
		dim *= sp->dim;
		sp = sp->next;
	}
        return dim;
}

int SizeProdSpace(ProdSpace h)
{
	int size = 0;
	for (; h != 0; h = h->next) {
		++size;
	}
	return size;
}
