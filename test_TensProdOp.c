#include "TensProdOp.h"
#include "TensProdOpImpl.h"

int test_embed()
{
  int errs = 0;
  TOp op = embed(0, 0, 0);
  return errs;
}

int main()
{
	int errs = 0;
  errs += test_embed();
	return errs;
}

