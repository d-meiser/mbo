#ifndef TENS_PROD_OP_H
#define TENS_PROD_OP_H

struct TOp;
struct ElemOp;

struct TOp *embed(const struct ElemOp *, const struct TSpace *);

#endif
