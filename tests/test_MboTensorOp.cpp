#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboAmplitude.h>
#include <vector>

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
  MboTensorOp A, B, C;
  struct MboAmplitude a, b, one, result, *arr, expectedResult;
  std::vector<struct MboAmplitude> x, y;
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

  /* x <- I * x + 0 * x */
  x.resize(mboProdSpaceDim(h2));
  std::fill(x.begin(), x.end(), one);
  y.resize(mboProdSpaceDim(h2));
  std::fill(y.begin(), y.end(), one);
  mboTensorOpIdentity(h2, &A);
  err = mboTensorOpMatVec(a, A, &x[0], b, &y[0], 0, mboProdSpaceDim(h2));
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(y[i].re, 1.0);
    EXPECT_FLOAT_EQ(y[i].im, 0.0);
  }
  mboTensorOpDestroy(&A);

  /* y <- a * I * x + b * y
   * expected result:
   * a * one * one + b * b */
  x.resize(mboProdSpaceDim(h2));
  a.re = 2.5;
  a.im = 22.0;
  std::fill(x.begin(), x.end(), one);
  b.re = 3.0;
  b.im = -1.7;
  y.resize(mboProdSpaceDim(h2));
  std::fill(y.begin(), y.end(), b);
  result.re = a.re + b.re * b.re - b.im * b.im;
  result.im = a.im + b.re * b.im + b.im * b.re;
  mboTensorOpIdentity(h2, &A);
  err = mboTensorOpMatVec(a, A, &x[0], b, &y[0], 0, mboProdSpaceDim(h2));
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(y[i].re, result.re);
    EXPECT_FLOAT_EQ(y[i].im, result.im);
  }
  mboTensorOpDestroy(&A);

  /* y <-  s_minus(0) * x */
  x.resize(mboProdSpaceDim(h2));
  a.re = 2.5;
  a.im = 22.0;
  std::fill(x.begin(), x.end(), a);
  b.re = 3.0;
  b.im = -1.7;
  y.resize(mboProdSpaceDim(h2));
  std::fill(y.begin(), y.end(), b);
  mboTensorOpNull(h2, &A);
  eop = mboSigmaMinus();
  mboTensorOpAddTo(eop, 0, A);
  b.re = 0.0;
  b.im = 0.0;
  err = mboTensorOpMatVec(one, A, &x[0], b, &y[0], 0, mboProdSpaceDim(h2));
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_FLOAT_EQ(a.re, y[i].re);
    EXPECT_FLOAT_EQ(a.im, y[i].im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(0, y[i].re);
    EXPECT_FLOAT_EQ(0, y[i].im);
  }
  mboTensorOpDestroy(&A);
  mboElemOpDestroy(&eop);

  /* y <-  (I + s_minus(0)) * x */
  x.resize(mboProdSpaceDim(h2));
  a.re = 2.5;
  a.im = 22.0;
  std::fill(x.begin(), x.end(), a);
  b.re = 3.0;
  b.im = -1.7;
  y.resize(mboProdSpaceDim(h2));
  std::fill(y.begin(), y.end(), b);
  mboTensorOpIdentity(h2, &A);
  eop = mboSigmaMinus();
  mboTensorOpAddTo(eop, 0, A);
  b.re = 0.0;
  b.im = 0.0;
  err = mboTensorOpMatVec(one, A, &x[0], b, &y[0], 0, mboProdSpaceDim(h2));
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_FLOAT_EQ(y[i].re, 2.0 * a.re);
    EXPECT_FLOAT_EQ(y[i].im, 2.0 * a.im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_FLOAT_EQ(y[i].re, a.re);
    EXPECT_FLOAT_EQ(y[i].im, a.im);
  }
  mboTensorOpDestroy(&A);
  mboElemOpDestroy(&eop);

  x.resize(mboProdSpaceDim(h2));
  a.re = 2.5;
  a.im = 22.0;
  std::fill(x.begin(), x.end(), a);
  b.re = 3.0;
  b.im = -1.7;
  y.resize(mboProdSpaceDim(h2));
  std::fill(y.begin(), y.end(), b);
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
  err = mboTensorOpMatVec(one, C, &x[0], b, &y[0], 0, mboProdSpaceDim(h2));
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
    EXPECT_FLOAT_EQ(y[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(y[i].im, expectedResult.im);
  }
  free(dims);
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

TEST(MboTensorOp, RowOffsetsNull) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(null, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(0, i[2]);

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsIdentity) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(id, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(1, i[1]);
  EXPECT_EQ(2, i[2]);

  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsIdentityEmptyRange) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);
  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(id, 3, 2, &i[0]);
  mboTensorOpRowOffsets(id, 2, 2, &i[0]);
  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsIdentitySubrange) {
  int dim = 5;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);

  std::vector<int> i(3);
  mboTensorOpRowOffsets(id, 2, 4, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(1, i[1]);
  EXPECT_EQ(2, i[2]);

  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(Sp, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(1, i[2]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsInBiggerSpace) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Sp);

  MboGlobInd totalDim = mboProdSpaceDim(h);
  std::vector<int> i(totalDim + 1);
  mboTensorOpRowOffsets(Sp, 0, totalDim, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(0, i[2]);
  EXPECT_EQ(0, i[3]);
  EXPECT_EQ(0, i[4]);
  EXPECT_EQ(1, i[5]);
  EXPECT_EQ(2, i[6]);
  EXPECT_EQ(3, i[7]);
  EXPECT_EQ(4, i[8]);
  EXPECT_EQ(4, i[9]);
  EXPECT_EQ(4, i[10]);
  EXPECT_EQ(4, i[11]);
  EXPECT_EQ(4, i[12]);
  EXPECT_EQ(5, i[13]);
  EXPECT_EQ(6, i[14]);
  EXPECT_EQ(7, i[15]);
  EXPECT_EQ(8, i[16]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsTwoEmbeddings) {
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
  std::vector<int> i(totalDim + 1);
  mboTensorOpRowOffsets(Sp, 0, totalDim, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(1, i[2]);
  EXPECT_EQ(1, i[3]);
  EXPECT_EQ(2, i[4]);
  EXPECT_EQ(3, i[5]);
  EXPECT_EQ(5, i[6]);
  EXPECT_EQ(6, i[7]);
  EXPECT_EQ(8, i[8]);
  EXPECT_EQ(8, i[9]);
  EXPECT_EQ(9, i[10]);
  EXPECT_EQ(9, i[11]);
  EXPECT_EQ(10, i[12]);
  EXPECT_EQ(11, i[13]);
  EXPECT_EQ(13, i[14]);
  EXPECT_EQ(14, i[15]);
  EXPECT_EQ(16, i[16]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsSubRange) {
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
  std::vector<int> i(4);
  mboTensorOpRowOffsets(Sp, 5, 8, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(2, i[1]);
  EXPECT_EQ(3, i[2]);
  EXPECT_EQ(5, i[3]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, SparseMatrixNull) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(null, 0, 2, &i[0]);
  std::vector<int> j(i[dim]);
  std::vector<struct MboAmplitude> a(i[dim]);
  mboTensorOpSparseMatrix(null, 0, 2, &i[0], &j[0], &a[0]);

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, SparseMatrixIdentity) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(identity, 0, 2, &i[0]);
  std::vector<int> j(i[dim]);
  std::vector<struct MboAmplitude> a(i[dim]);
  mboTensorOpSparseMatrix(identity, 0, 2, &i[0], &j[0], &a[0]);
  EXPECT_EQ(0, j[0]);
  EXPECT_EQ(1, j[1]);
  EXPECT_FLOAT_EQ(1.0, a[0].re);
  EXPECT_FLOAT_EQ(0.0, a[0].im);
  EXPECT_FLOAT_EQ(1.0, a[1].re);
  EXPECT_FLOAT_EQ(0.0, a[1].im);

  mboTensorOpDestroy(&identity);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, SparseMatrixSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<int> i(dim + 1);
  mboTensorOpRowOffsets(Sp, 0, 2, &i[0]);
  std::vector<int> j(i[dim]);
  std::vector<struct MboAmplitude> a(i[dim]);
  mboTensorOpSparseMatrix(Sp, 0, 2, &i[0], &j[0], &a[0]);
  EXPECT_EQ(0, j[0]);
  EXPECT_FLOAT_EQ(1.0, a[0].re);
  EXPECT_FLOAT_EQ(0.0, a[0].im);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, SparseMatrixInBiggerSpace) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Sp);

  MboGlobInd totalDim = mboProdSpaceDim(h);
  std::vector<int> i(totalDim + 1);
  mboTensorOpRowOffsets(Sp, 0, totalDim, &i[0]);
  std::vector<int> j(i[mboProdSpaceDim(h)]);
  std::vector<struct MboAmplitude> a(i[mboProdSpaceDim(h)]);
  mboTensorOpSparseMatrix(Sp, 0, mboProdSpaceDim(h), &i[0], &j[0], &a[0]);
  EXPECT_EQ(8, j.size());
  EXPECT_EQ(0, j[0]);
  EXPECT_EQ(1, j[1]);
  EXPECT_EQ(2, j[2]);
  EXPECT_EQ(3, j[3]);
  EXPECT_EQ(8, j[4]);
  EXPECT_EQ(9, j[5]);
  EXPECT_EQ(10, j[6]);
  EXPECT_EQ(11, j[7]);
  EXPECT_EQ(8, a.size());
  for (size_t i = 0; i < a.size(); ++i) {
	  EXPECT_FLOAT_EQ(1.0, a[i].re) << "i = " << i;
	  EXPECT_FLOAT_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

static void printMatrixPattern(MboTensorOp op) {
  MboGlobInd dim = mboProdSpaceDim(mboTensorOpGetSpace(op));
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboTensorOpDenseMatrix(op, &mat[0]);
  std::cout << "\n\n";
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      if (mat[i * dim + j].re * mat[i * dim + j].re +
          mat[i * dim + j].im * mat[i * dim + j].im > 1.0e-16) {
	      std::cout << "x";
      } else {
	      std::cout << ".";
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n\n";
}

TEST(MboTensorOp, SparseMatrixTwoEmbeddings) {
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
  std::vector<int> i(totalDim + 1);
  mboTensorOpRowOffsets(Sp, 0, totalDim, &i[0]);
  std::vector<int> j(i[mboProdSpaceDim(h)]);
  std::vector<struct MboAmplitude> a(i[mboProdSpaceDim(h)]);
  mboTensorOpSparseMatrix(Sp, 0, mboProdSpaceDim(h), &i[0], &j[0], &a[0]);
  EXPECT_EQ(16, j.size());
  EXPECT_EQ(0, j[0]);
  EXPECT_EQ(2, j[1]);
  EXPECT_EQ(0, j[2]);
  EXPECT_EQ(1, j[3]);
  EXPECT_EQ(4, j[4]);
  EXPECT_EQ(2, j[5]);
  EXPECT_EQ(3, j[6]);
  EXPECT_EQ(6, j[7]);
  EXPECT_EQ(8, j[8]);
  EXPECT_EQ(10, j[9]);
  EXPECT_EQ(8, j[10]);
  EXPECT_EQ(9, j[11]);
  EXPECT_EQ(12, j[12]);
  EXPECT_EQ(10, j[13]);
  EXPECT_EQ(11, j[14]);
  EXPECT_EQ(14, j[15]);
  EXPECT_EQ(16, a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_FLOAT_EQ(1.0, a[i].re) << "i = " << i;
    EXPECT_FLOAT_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, SparseMatrixTwoEmbeddingsSubrange) {
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
  std::vector<int> i(totalDim + 1);
  MboGlobInd rmin = 3;
  MboGlobInd rmax = 11;
  mboTensorOpRowOffsets(Sp, rmin, rmax, &i[0]);
  std::vector<int> j(i[rmax - rmin]);
  std::vector<struct MboAmplitude> a(i[rmax - rmin]);
  mboTensorOpSparseMatrix(Sp, rmin, rmax, &i[0], &j[0], &a[0]);
  EXPECT_EQ(8, j.size());
  EXPECT_EQ(2, j[0]);
  EXPECT_EQ(0, j[1]);
  EXPECT_EQ(1, j[2]);
  EXPECT_EQ(4, j[3]);
  EXPECT_EQ(2, j[4]);
  EXPECT_EQ(3, j[5]);
  EXPECT_EQ(6, j[6]);
  EXPECT_EQ(8, j[7]);
  EXPECT_EQ(8, a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_FLOAT_EQ(1.0, a[i].re) << "i = " << i;
    EXPECT_FLOAT_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetDiagonalNull) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<struct MboAmplitude> diag(dim);
  mboTensorOpDiagonal(null, 0, dim, &diag[0]);

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetDiagonalIdentity) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  std::vector<struct MboAmplitude> diag(dim);
  mboTensorOpDiagonal(identity, 0, dim, &diag[0]);
  for (int i = 0; i < dim; ++i) {
    EXPECT_FLOAT_EQ(1.0, diag[i].re) << "i = " << i;
    EXPECT_FLOAT_EQ(0.0, diag[i].im) << "i = " << i;
  }

  mboTensorOpDestroy(&identity);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetDiagonalSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<struct MboAmplitude> diag(dim);
  mboTensorOpDiagonal(Sp, 0, dim, &diag[0]);
  EXPECT_FLOAT_EQ(0.0, diag[0].re);
  EXPECT_FLOAT_EQ(0.0, diag[0].im);
  EXPECT_FLOAT_EQ(0.0, diag[1].re);
  EXPECT_FLOAT_EQ(0.0, diag[1].im);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetDiagonalSigmaZ) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sz;
  mboTensorOpNull(h, &Sz);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 0, Sz);

  std::vector<struct MboAmplitude> diag(dim);
  mboTensorOpDiagonal(Sz, 0, dim, &diag[0]);
  EXPECT_FLOAT_EQ(-1.0, diag[0].re);
  EXPECT_FLOAT_EQ(0.0, diag[0].im);
  EXPECT_FLOAT_EQ(1.0, diag[1].re);
  EXPECT_FLOAT_EQ(0.0, diag[1].im);

  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&Sz);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, GetDiagonalSigmaZInBiggerSpace) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sz;
  mboTensorOpNull(h, &Sz);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 2, Sz);

  std::vector<struct MboAmplitude> diag(mboProdSpaceDim(h));
  mboTensorOpDiagonal(Sz, 0, mboProdSpaceDim(h), &diag[0]);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (j / 2 == 0) {
        EXPECT_FLOAT_EQ(-1.0, diag[i * 4 + j].re) << " i == " << i
                                                  << " j == " << j;
        EXPECT_FLOAT_EQ(0.0, diag[i * 4 + j].im) << " i == " << i
                                                 << " j == " << j;
      } else {
        EXPECT_FLOAT_EQ(1.0, diag[i * 4 + j].re) << " i == " << i
                                                 << " j == " << j;
        EXPECT_FLOAT_EQ(0.0, diag[i * 4 + j].im) << " i == " << i
                                                 << " j == " << j;
      }
    }
  }

  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&Sz);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DeleteDiagonalNull) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  mboTensorOpDeleteDiagonal(null);
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboTensorOpDenseMatrix(null, &mat[0]);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_FLOAT_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_FLOAT_EQ(0, mat[i].im) << " i == " << i;
  }

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DeleteDiagonalIdentity) {
  int dim = 3;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  mboTensorOpDeleteDiagonal(identity);
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboTensorOpDenseMatrix(identity, &mat[0]);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_FLOAT_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_FLOAT_EQ(0, mat[i].im) << " i == " << i;
  }

  mboTensorOpDestroy(&identity);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DeleteDiagonalSz) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Op;
  mboTensorOpNull(h, &Op);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 1, Op);

  mboTensorOpDeleteDiagonal(Op);
  int D = mboProdSpaceDim(h);
  std::vector<struct MboAmplitude> mat(D * D);
  mboTensorOpDenseMatrix(Op, &mat[0]);
  for (int i = 0; i < D * D; ++i) {
    EXPECT_FLOAT_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_FLOAT_EQ(0, mat[i].im) << " i == " << i;
  }

  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&Op);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DeleteDiagonalSp) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Op;
  mboTensorOpNull(h, &Op);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Op);

  mboTensorOpDeleteDiagonal(Op);
  int D = mboProdSpaceDim(h);
  std::vector<int> I(D + 1);
  mboTensorOpRowOffsets(Op, 0, D, &I[0]);
  EXPECT_EQ(0, I[0]);
  EXPECT_EQ(0, I[1]);
  EXPECT_EQ(0, I[2]);
  EXPECT_EQ(0, I[3]);
  EXPECT_EQ(0, I[4]);
  EXPECT_EQ(1, I[5]);
  EXPECT_EQ(2, I[6]);
  EXPECT_EQ(3, I[7]);
  EXPECT_EQ(4, I[8]);
  EXPECT_EQ(4, I[9]);
  EXPECT_EQ(4, I[10]);
  EXPECT_EQ(4, I[11]);
  EXPECT_EQ(4, I[12]);
  EXPECT_EQ(5, I[13]);
  EXPECT_EQ(6, I[14]);
  EXPECT_EQ(7, I[15]);
  EXPECT_EQ(8, I[16]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Op);
  mboProdSpaceDestroy(&h);
}
