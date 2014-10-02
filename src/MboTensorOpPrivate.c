#include <MboTensorOpPrivate.h>

int mboTensorOpGetNumTerms(MboTensorOp op)
{
  return op->numTerms;
}

struct SimpleTOp* mboTensorOpGetSimpleTOps(MboTensorOp op)
{
  return op->sum;
}

