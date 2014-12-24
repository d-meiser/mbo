#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
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
  EXPECT_NE(0, mboProdSpaceEqual(h, h2));
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

TEST(MboTensorOp, Copy) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboElemOp eop = mboSigmaZ();
  MboTensorOp A;
  mboTensorOpNull(h, &A);
  mboTensorOpAddTo(eop, 0, A);

  MboTensorOp B = mboTensorOpCopy(A);
  MboNumOp Acomp = mboNumOpCompile(A);
  MboNumOp Bcomp = mboNumOpCompile(B);
  std::vector<struct MboAmplitude> x(2);
  x[0].re = 134.3;
  x[0].im = 11.0;
  x[1].re = -12.0;
  x[1].im = 0.57;
  struct MboAmplitude one;
  one.re = 1.0;
  one.im = 0.0;
  struct MboAmplitude zero;
  zero.re = 0.0;
  zero.im = 0.0;
  std::vector<struct MboAmplitude> ya(2);
  std::vector<struct MboAmplitude> yb(2);
  mboNumOpMatVec(one, Acomp, &x[0], zero, &ya[0]);
  mboNumOpMatVec(one, Bcomp, &x[0], zero, &yb[0]);
  EXPECT_DOUBLE_EQ(ya[0].re, yb[0].re);
  EXPECT_DOUBLE_EQ(ya[0].im, yb[0].im);
  EXPECT_DOUBLE_EQ(ya[1].re, yb[1].re);
  EXPECT_DOUBLE_EQ(ya[1].im, yb[1].im);

  mboNumOpDestroy(&Acomp);
  mboNumOpDestroy(&Bcomp);
  mboTensorOpDestroy(&A);
  mboTensorOpDestroy(&B);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h);
}

