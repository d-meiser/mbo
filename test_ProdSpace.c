#include "ProdSpace.h"

int testCreateProdSpace()
{
  int errs = 0;
  ProdSpace sp;
  sp = CreateProdSpace(2);
  DestroyProdSpace(sp);
  return errs;
}

int main()
{
	int errs = 0;
        errs += testCreateProdSpace();
	return errs;
}
