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
