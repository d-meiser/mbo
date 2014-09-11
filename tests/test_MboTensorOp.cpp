#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboAmplitude.h>

TEST(MboTensorOp, Null) {
  MboTensorOp op;
  MboProdSpace h = mboProdSpaceCreate(1);
  mboTensorOpNull(h, &op);
  EXPECT_NE(op, (MboTensorOp)0);
  EXPECT_EQ(mboTensorOpCheck(op), 0);
  mboTensorOpDestroy(&op);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, Identity) {
  MboTensorOp op;
  MboProdSpace h = mboProdSpaceCreate(1);
  mboTensorOpIdentity(h, &op);
  EXPECT_NE(op, (MboTensorOp)0);
  EXPECT_EQ(mboTensorOpCheck(op), 0);
  mboTensorOpDestroy(&op);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetSpace) {
  MboTensorOp op;
  MboProdSpace h = mboProdSpaceCreate(5);
  mboTensorOpIdentity(h, &op);
  MboProdSpace h2 = mboTensorOpGetSpace(op);
  EXPECT_TRUE(mboProdSpaceEqual(h, h2));
  mboTensorOpDestroy(&op);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, AddTo) {
  MboTensorOp op;
  MboElemOp eop;
  MboProdSpace h;
  struct MboAmplitude alpha;

  h = mboProdSpaceCreate(20);

  mboElemOpCreate(&eop);
  alpha.re = 5.4;
  alpha.im = 2.9;
  mboElemOpAddTo(14, 15, &alpha, &eop);

  mboTensorOpNull(h, &op);
  mboTensorOpAddTo(eop, 0, op);
  EXPECT_EQ(mboTensorOpCheck(op), 0);

  mboTensorOpDestroy(&op);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, AddScaledTo) {
  MboTensorOp op;
  MboElemOp eop;
  MboProdSpace h;
  struct MboAmplitude alpha;

  h = mboProdSpaceCreate(20);

  mboElemOpCreate(&eop);
  alpha.re = 5.4;
  alpha.im = 2.9;
  mboElemOpAddTo(14, 15, &alpha, &eop);

  mboTensorOpNull(h, &op);
  mboTensorOpAddScaledTo(&alpha, eop, 0, op);
  EXPECT_EQ(mboTensorOpCheck(op), 0);

  mboTensorOpDestroy(&op);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, Mul) {
  int i;
  int N = 10;
  MboTensorOp op1, op2, op3;
  MboElemOp eop1, eop2;
  MboProdSpace h1, h2;
  struct MboAmplitude alpha, beta;

  alpha.re = 20;
  alpha.im = 30;
  beta.re = 15;
  beta.im = -19;

  h1 = mboProdSpaceCreate(5);
  h2 = mboProdSpaceCreate(0);
  for (i = 0; i < N; ++i) {
    mboProdSpaceMul(h1, &h2);
  }

  mboElemOpCreate(&eop1);
  mboElemOpAddTo(14, 15, &alpha, &eop1);
  mboElemOpAddTo(1, 5, &beta, &eop1);
  mboElemOpAddTo(2, 3, &beta, &eop1);

  mboElemOpCreate(&eop2);
  mboElemOpAddTo(8, 3, &beta, &eop2);
  mboElemOpAddTo(3, 4, &alpha, &eop2);

  mboTensorOpNull(h2, &op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboTensorOpNull(h2, &op1);
  mboTensorOpAddTo(eop1, 0, op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboTensorOpNull(h2, &op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpAddTo(eop1, 0, op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboTensorOpNull(h2, &op1);
  mboTensorOpAddTo(eop1, 0, op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpAddTo(eop2, 0, op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboTensorOpNull(h2, &op1);
  mboTensorOpAddTo(eop1, 0, op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpAddTo(eop2, 1, op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboTensorOpNull(h2, &op1);
  mboTensorOpAddTo(eop1, 0, op1);
  mboTensorOpAddTo(eop1, 1, op1);
  mboTensorOpNull(h2, &op2);
  mboTensorOpAddTo(eop2, 1, op2);
  mboTensorOpAddTo(eop2, 2, op2);
  mboTensorOpNull(h2, &op3);
  mboTensorOpMul(op1, op2, &op3);
  EXPECT_EQ(mboTensorOpCheck(op3), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);
  mboTensorOpDestroy(&op3);

  mboElemOpDestroy(&eop1);
  mboElemOpDestroy(&eop2);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);
}

TEST(MboTensorOp, Plus) {
  int N = 5;
  MboProdSpace hTot, h1;
  MboTensorOp op1, op2;
  MboElemOp sz, sp, sm;

  hTot = mboProdSpaceCreate(0);
  h1 = mboProdSpaceCreate(2);
  for (int i = 0; i < N; ++i) {
    mboProdSpaceMul(h1, &hTot);
  }

  sz = mboSigmaZ();
  sp = mboSigmaPlus();
  sm = mboSigmaMinus();

  mboTensorOpNull(hTot, &op1);
  mboTensorOpNull(hTot, &op2);
  mboTensorOpPlus(op1, &op2);
  EXPECT_EQ(mboTensorOpCheck(op2), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);

  mboTensorOpNull(hTot, &op1);
  mboTensorOpNull(hTot, &op2);
  mboTensorOpAddTo(sz, 1, op2);
  mboTensorOpPlus(op1, &op2);
  EXPECT_EQ(mboTensorOpCheck(op2), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);

  mboTensorOpNull(hTot, &op1);
  mboTensorOpAddTo(sz, 1, op1);
  mboTensorOpNull(hTot, &op2);
  mboTensorOpPlus(op1, &op2);
  EXPECT_EQ(mboTensorOpCheck(op2), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);

  mboTensorOpNull(hTot, &op1);
  mboTensorOpNull(hTot, &op2);
  mboTensorOpAddTo(sz, 0, op1);
  mboTensorOpAddTo(sz, 1, op1);
  mboTensorOpAddTo(sp, 0, op2);
  mboTensorOpAddTo(sp, 2, op2);
  mboTensorOpAddTo(sp, 3, op2);
  mboTensorOpAddTo(sm, 3, op2);
  mboTensorOpAddTo(sm, 2, op2);
  EXPECT_EQ(mboTensorOpCheck(op2), 0);
  mboTensorOpDestroy(&op1);
  mboTensorOpDestroy(&op2);

  mboElemOpDestroy(&sz);
  mboElemOpDestroy(&sp);
  mboElemOpDestroy(&sm);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);
}

TEST(MboTensorOp, Scale) {
  MboTensorOp a;
  MboProdSpace h;
  struct MboAmplitude alpha;

  alpha.re = 2.0;
  alpha.im = -55.55;

  h = mboProdSpaceCreate(0);
  mboTensorOpIdentity(h, &a);
  mboTensorOpScale(&alpha, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);

  h = mboProdSpaceCreate(2);
  mboTensorOpNull(h, &a);
  mboTensorOpScale(&alpha, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);

  h = mboProdSpaceCreate(2);
  mboTensorOpIdentity(h, &a);
  mboTensorOpScale(&alpha, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, Kron) {
  MboTensorOp *ops, result;
  MboProdSpace h1, h2, h3, hTot;
  MBO_STATUS err;
  MboElemOp sz;
  MboLocInd dims[3];

  h1 = mboProdSpaceCreate(2);
  hTot = mboProdSpaceCreate(0);
  mboTensorOpNull(hTot, &result);
  int n = 1;
  ops = (MboTensorOp *)malloc(n * sizeof(*ops));
  mboTensorOpIdentity(h1, &ops[0]);
  err = mboTensorOpKron(1, ops, &result);
  EXPECT_EQ(err, MBO_SPACE_MISMATCH);
  mboTensorOpDestroy(&result);
  mboTensorOpDestroy(&ops[0]);
  free(ops);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);

  h1 = mboProdSpaceCreate(2);
  hTot = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &hTot);
  mboTensorOpNull(hTot, &result);
  n = 1;
  ops = (MboTensorOp *)malloc(n * sizeof(*ops));
  mboTensorOpIdentity(h1, &ops[0]);
  err = mboTensorOpKron(n, ops, &result);
  EXPECT_EQ(err, MBO_SUCCESS);
  mboTensorOpDestroy(&ops[0]);
  free(ops);
  mboTensorOpDestroy(&result);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);

  h1 = mboProdSpaceCreate(2);
  h2 = mboProdSpaceCreate(3);
  hTot = mboProdSpaceCreate(0);
  mboProdSpaceMul(h2, &hTot);
  mboProdSpaceMul(h1, &hTot);
  mboTensorOpNull(hTot, &result);
  n = 2;
  ops = (MboTensorOp *)malloc(n * sizeof(*ops));
  mboTensorOpIdentity(h1, &ops[0]);
  mboTensorOpIdentity(h2, &ops[1]);
  err = mboTensorOpKron(n, ops, &result);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_EQ(mboTensorOpCheck(result), MBO_SUCCESS);
  mboTensorOpDestroy(&ops[0]);
  mboTensorOpDestroy(&ops[1]);
  free(ops);
  mboTensorOpDestroy(&result);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);

  sz = mboSigmaZ();
  h1 = mboProdSpaceCreate(2);
  h2 = mboProdSpaceCreate(3);
  hTot = mboProdSpaceCreate(0);
  mboProdSpaceMul(h2, &hTot);
  mboProdSpaceMul(h1, &hTot);
  mboTensorOpNull(hTot, &result);
  n = 2;
  ops = (MboTensorOp *)malloc(n * sizeof(*ops));
  mboTensorOpIdentity(h1, &ops[0]);
  mboTensorOpIdentity(h2, &ops[1]);
  mboTensorOpAddTo(sz, 0, ops[0]);
  err = mboTensorOpKron(n, ops, &result);
  EXPECT_EQ(err, MBO_SUCCESS);
  mboTensorOpDestroy(&ops[0]);
  mboTensorOpDestroy(&ops[1]);
  free(ops);
  mboTensorOpDestroy(&result);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);
  mboElemOpDestroy(&sz);

  sz = mboSigmaZ();
  h1 = mboProdSpaceCreate(2);
  h2 = mboProdSpaceCreate(3);
  h3 = mboProdSpaceCreate(2);
  hTot = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &hTot);
  mboProdSpaceMul(h2, &hTot);
  mboProdSpaceMul(h3, &hTot);
  mboTensorOpNull(hTot, &result);
  n = 3;
  ops = (MboTensorOp *)malloc(n * sizeof(*ops));
  mboTensorOpIdentity(h1, &ops[0]);
  mboTensorOpIdentity(h2, &ops[1]);
  mboTensorOpIdentity(h3, &ops[2]);
  mboTensorOpAddTo(sz, 0, ops[0]);
  mboTensorOpAddTo(sz, 0, ops[1]);
  mboTensorOpAddTo(sz, 0, ops[2]);
  err = mboTensorOpKron(n, ops, &result);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_EQ(mboTensorOpCheck(result), 0);
  mboTensorOpDestroy(&ops[0]);
  mboTensorOpDestroy(&ops[1]);
  mboTensorOpDestroy(&ops[2]);
  free(ops);
  mboTensorOpDestroy(&result);
  mboProdSpaceDestroy(&hTot);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);
  mboProdSpaceDestroy(&h3);
  mboElemOpDestroy(&sz);
}

static MboGlobInd computeBlockSize(int N, MboLocInd *dims) {
  int i;
  MboGlobInd blockSize = 1;
  for (i = 0; i < N; ++i) {
    blockSize *= (MboGlobInd)dims[i];
  }
  return blockSize;
}

TEST(MboTensorOp, MatVec) {
  MboGlobInd i;
  MboLocInd *dims;
  MboProdSpace h1, h2;
  MboVec x, y;
  MboTensorOp A, B, C;
  struct MboAmplitude a, b, one, result, *arr, expectedResult;
  MboElemOp eop;
  MBO_STATUS err;

  one.re = 1.0;
  one.im = 0.0;
  a = one;
  b.re = 0.0;
  b.im = 0.0;

  /* set up spaces */
  h1 = mboProdSpaceCreate(2);
  h2 = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceDestroy(&h1);
  h1 = mboProdSpaceCreate(3);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceDestroy(&h1);
  h1 = mboProdSpaceCreate(2);
  mboProdSpaceMul(h1, &h2);

  /* mismatching dimensions */
  mboVecCreate(1l + mboProdSpaceDim(h2), &x);
  mboTensorOpIdentity(h2, &A);
  err = mboTensorOpMatVec(&a, A, x, &b, x);
  EXPECT_EQ(err, MBO_DIMENSIONS_MISMATCH);
  mboTensorOpDestroy(&A);
  mboVecDestroy(&x);

  /* x <- I * x + 0 * x */
  mboVecCreate(mboProdSpaceDim(h2), &x);
  mboVecSet(&one, x);
  mboTensorOpIdentity(h2, &A);
  err = mboTensorOpMatVec(&a, A, x, &b, x);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(x, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, 1.0);
    EXPECT_FLOAT_EQ(arr[i].im, 0.0);
  }
  mboTensorOpDestroy(&A);
  mboVecDestroy(&x);

  /* y <- a * I * x + b * y
   * expected result:
   * a * one * one + b * b */
  mboVecCreate(mboProdSpaceDim(h2), &x);
  a.re = 2.5;
  a.im = 22.0;
  mboVecSet(&one, x);
  mboVecCreate(mboProdSpaceDim(h2), &y);
  b.re = 3.0;
  b.im = -1.7;
  mboVecSet(&b, y);
  result.re = a.re + b.re * b.re - b.im * b.im;
  result.im = a.im + b.re * b.im + b.im * b.re;
  mboTensorOpIdentity(h2, &A);
  err = mboTensorOpMatVec(&a, A, x, &b, y);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(y, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, result.re);
    EXPECT_FLOAT_EQ(arr[i].im, result.im);
  }
  mboTensorOpDestroy(&A);
  mboVecDestroy(&x);
  mboVecDestroy(&y);

  /* y <-  s_minus(0) * x */
  mboVecCreate(mboProdSpaceDim(h2), &x);
  a.re = 2.5;
  a.im = 22.0;
  mboVecSet(&a, x);
  mboVecCreate(mboProdSpaceDim(h2), &y);
  b.re = 3.0;
  b.im = -1.7;
  mboVecSet(&b, y);
  mboTensorOpNull(h2, &A);
  eop = mboSigmaMinus();
  mboTensorOpAddTo(eop, 0, A);
  b.re = 0.0;
  b.im = 0.0;
  err = mboTensorOpMatVec(&one, A, x, &b, y);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(y, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, a.re);
    EXPECT_FLOAT_EQ(arr[i].im, a.im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, 0);
    EXPECT_FLOAT_EQ(arr[i].im, 0);
  }
  mboTensorOpDestroy(&A);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
  mboElemOpDestroy(&eop);

  /* y <-  (I + s_minus(0)) * x */
  mboVecCreate(mboProdSpaceDim(h2), &x);
  a.re = 2.5;
  a.im = 22.0;
  mboVecSet(&a, x);
  mboVecCreate(mboProdSpaceDim(h2), &y);
  b.re = 3.0;
  b.im = -1.7;
  mboVecSet(&b, y);
  mboTensorOpIdentity(h2, &A);
  eop = mboSigmaMinus();
  mboTensorOpAddTo(eop, 0, A);
  b.re = 0.0;
  b.im = 0.0;
  err = mboTensorOpMatVec(&one, A, x, &b, y);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(y, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, 2.0 * a.re);
    EXPECT_FLOAT_EQ(arr[i].im, 2.0 * a.im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, a.re);
    EXPECT_FLOAT_EQ(arr[i].im, a.im);
  }
  mboTensorOpDestroy(&A);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
  mboElemOpDestroy(&eop);

  mboVecCreate(mboProdSpaceDim(h2), &x);
  a.re = 2.5;
  a.im = 22.0;
  mboVecSet(&a, x);
  mboVecCreate(mboProdSpaceDim(h2), &y);
  b.re = 3.0;
  b.im = -1.7;
  mboVecSet(&b, y);
  mboTensorOpNull(h2, &A);
  eop = mboSigmaMinus();
  mboTensorOpAddTo(eop, 0, A);
  mboElemOpDestroy(&eop);
  mboTensorOpNull(h2, &B);
  eop = mboSigmaPlus();
  mboTensorOpAddTo(eop, 1, B);
  mboElemOpDestroy(&eop);
  mboTensorOpNull(h2, &C);
  mboTensorOpMul(A, B, &C);
  mboTensorOpDestroy(&A);
  mboTensorOpDestroy(&B);
  err = mboTensorOpMatVec(&one, C, x, &b, y);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(y, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  dims = (MboLocInd *)malloc(3 * sizeof(*dims));
  mboProdSpaceGetDims(h2, 3, dims);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    expectedResult.re = b.re * b.re - b.im * b.im;
    expectedResult.im = b.re * b.im + b.im * b.re;
    if (((i / computeBlockSize(2, dims + 1)) % dims[0] == 0) &&
        ((i / computeBlockSize(1, dims + 2)) % dims[1] == 1)) {
      expectedResult.re += a.re;
      expectedResult.im += a.im;
    }
    EXPECT_FLOAT_EQ(arr[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(arr[i].im, expectedResult.im);
  }
  free(dims);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
  mboTensorOpDestroy(&C);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);
}

TEST(MboTensorOp, Flops) {
  MboTensorOp a;
  MboElemOp sz;
  MboProdSpace h, h1;
  double flops;

  h = mboProdSpaceCreate(2);
  mboTensorOpNull(h, &a);
  flops = mboTensorOpFlops(a);
  EXPECT_FLOAT_EQ(flops, 0.0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);

  h = mboProdSpaceCreate(5);
  mboTensorOpIdentity(h, &a);
  flops = mboTensorOpFlops(a);
  EXPECT_FLOAT_EQ(flops, 5 * 8);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);

  h1 = mboProdSpaceCreate(5);
  h = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboTensorOpIdentity(h, &a);
  flops = mboTensorOpFlops(a);
  EXPECT_FLOAT_EQ(flops, 5 * 5 * 5 * 8);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);
  mboProdSpaceDestroy(&h1);

  h1 = mboProdSpaceCreate(5);
  h = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboTensorOpIdentity(h, &a);
  sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 0, a);
  flops = mboTensorOpFlops(a);
  EXPECT_FLOAT_EQ(flops, (5 * 5 * 5 + 2 * 5 * 5) * 8);
  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);
  mboProdSpaceDestroy(&h1);
}

TEST(MboTensorOp, Check) {
  MboTensorOp a;
  MboProdSpace h;

  h = mboProdSpaceCreate(0);
  mboTensorOpNull(h, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboTensorOpIdentity(h, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);

  h = mboProdSpaceCreate(30);
  mboTensorOpNull(h, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboTensorOpIdentity(h, &a);
  EXPECT_EQ(mboTensorOpCheck(a), 0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DenseMatrixNull) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  MboGlobInd dim = mboProdSpaceDim(h);
  struct MboAmplitude *mat = new struct MboAmplitude[dim * dim];
  mboTensorOpDenseMatrix(null, mat);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_FLOAT_EQ(mat[i].re, 0);
    EXPECT_FLOAT_EQ(mat[i].im, 0);
  }

  delete [] mat;
  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DenseMatrixIdentity) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);

  MboGlobInd dim = mboProdSpaceDim(h);
  struct MboAmplitude *mat = new struct MboAmplitude[dim * dim];
  mboTensorOpDenseMatrix(id, mat);
  for (int i = 0; i < dim * dim; ++i) {
    struct MboAmplitude expectedResult;
    if (i % mboProdSpaceDim(h) == i / mboProdSpaceDim(h)) {
      expectedResult.re = 1.0;
      expectedResult.im = 0.0;
    } else {
      expectedResult.re = 0.0;
      expectedResult.im = 0.0;
    }
    EXPECT_FLOAT_EQ(mat[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(mat[i].im, expectedResult.im);
  }

  delete [] mat;
  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowNull) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<int> nnzs(dim);
  mboTensorOpGetNonZerosPerRow(null, 0, 2, &nnzs[0]);
  EXPECT_EQ(0, nnzs[0]);
  EXPECT_EQ(0, nnzs[1]);

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowIdentity) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);

  std::vector<int> nnzs(dim);
  mboTensorOpGetNonZerosPerRow(id, 0, 2, &nnzs[0]);
  EXPECT_EQ(1, nnzs[0]);
  EXPECT_EQ(1, nnzs[1]);

  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowIdentityEmptyRange) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);
  std::vector<int> nnzs(dim);
  mboTensorOpGetNonZerosPerRow(id, 3, 2, &nnzs[0]);
  mboTensorOpGetNonZerosPerRow(id, 2, 2, &nnzs[0]);
  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowIdentitySubrange) {
  int dim = 5;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);

  std::vector<int> nnzs(2);
  mboTensorOpGetNonZerosPerRow(id, 2, 4, &nnzs[0]);
  EXPECT_EQ(1, nnzs[0]);
  EXPECT_EQ(1, nnzs[1]);

  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<int> nnzs(dim);
  mboTensorOpGetNonZerosPerRow(Sp, 0, 2, &nnzs[0]);
  EXPECT_EQ(0, nnzs[0]);
  EXPECT_EQ(1, nnzs[1]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowInBiggerSpace) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Sp);

  MboGlobInd totalDim = mboProdSpaceDim(h);
  std::vector<int> nnzs(totalDim);
  mboTensorOpGetNonZerosPerRow(Sp, 0, totalDim, &nnzs[0]);
  EXPECT_EQ(0, nnzs[0]);
  EXPECT_EQ(0, nnzs[1]);
  EXPECT_EQ(0, nnzs[2]);
  EXPECT_EQ(0, nnzs[3]);
  EXPECT_EQ(1, nnzs[4]);
  EXPECT_EQ(1, nnzs[5]);
  EXPECT_EQ(1, nnzs[6]);
  EXPECT_EQ(1, nnzs[7]);
  EXPECT_EQ(0, nnzs[8]);
  EXPECT_EQ(0, nnzs[9]);
  EXPECT_EQ(0, nnzs[10]);
  EXPECT_EQ(0, nnzs[11]);
  EXPECT_EQ(1, nnzs[12]);
  EXPECT_EQ(1, nnzs[13]);
  EXPECT_EQ(1, nnzs[14]);
  EXPECT_EQ(1, nnzs[15]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowTwoEmbeddings) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Sp);
  mboTensorOpAddTo(sp, 3, Sp);

  MboGlobInd totalDim = mboProdSpaceDim(h);
  std::vector<int> nnzs(totalDim);
  mboTensorOpGetNonZerosPerRow(Sp, 0, totalDim, &nnzs[0]);
  EXPECT_EQ(0, nnzs[0]);
  EXPECT_EQ(1, nnzs[1]);
  EXPECT_EQ(0, nnzs[2]);
  EXPECT_EQ(1, nnzs[3]);
  EXPECT_EQ(1, nnzs[4]);
  EXPECT_EQ(2, nnzs[5]);
  EXPECT_EQ(1, nnzs[6]);
  EXPECT_EQ(2, nnzs[7]);
  EXPECT_EQ(0, nnzs[8]);
  EXPECT_EQ(1, nnzs[9]);
  EXPECT_EQ(0, nnzs[10]);
  EXPECT_EQ(1, nnzs[11]);
  EXPECT_EQ(1, nnzs[12]);
  EXPECT_EQ(2, nnzs[13]);
  EXPECT_EQ(1, nnzs[14]);
  EXPECT_EQ(2, nnzs[15]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetNonZerosPerRowSubRange) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Sp);
  mboTensorOpAddTo(sp, 3, Sp);

  MboGlobInd totalDim = mboProdSpaceDim(h);
  std::vector<int> nnzs(3);
  mboTensorOpGetNonZerosPerRow(Sp, 5, 8, &nnzs[0]);
  EXPECT_EQ(2, nnzs[0]);
  EXPECT_EQ(1, nnzs[1]);
  EXPECT_EQ(2, nnzs[2]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}
