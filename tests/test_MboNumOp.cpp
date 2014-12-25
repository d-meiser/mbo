#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
#include <MboAmplitude.h>
#include <vector>
#include <algorithm>

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
  MboNumOp Ac;
  struct MboAmplitude a, b, one, result, expectedResult;
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
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(A, &Ac));
  err = mboNumOpMatVec(a, Ac, &x[0], b, &y[0]);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_DOUBLE_EQ(y[i].re, 1.0);
    EXPECT_DOUBLE_EQ(y[i].im, 0.0);
  }
  mboTensorOpDestroy(&A);
  mboNumOpDestroy(&Ac);

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
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(A, &Ac));
  err = mboNumOpMatVec(a, Ac, &x[0], b, &y[0]);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_DOUBLE_EQ(y[i].re, result.re);
    EXPECT_DOUBLE_EQ(y[i].im, result.im);
  }
  mboTensorOpDestroy(&A);
  mboNumOpDestroy(&Ac);

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
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(A, &Ac));
  err = mboNumOpMatVec(one, Ac, &x[0], b, &y[0]);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_DOUBLE_EQ(a.re, y[i].re);
    EXPECT_DOUBLE_EQ(a.im, y[i].im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_DOUBLE_EQ(0, y[i].re);
    EXPECT_DOUBLE_EQ(0, y[i].im);
  }
  mboTensorOpDestroy(&A);
  mboNumOpDestroy(&Ac);
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
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(A, &Ac));
  err = mboNumOpMatVec(one, Ac, &x[0], b, &y[0]);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (i = 0; i < mboProdSpaceDim(h2) / 2; ++i) {
    EXPECT_DOUBLE_EQ(y[i].re, 2.0 * a.re);
    EXPECT_DOUBLE_EQ(y[i].im, 2.0 * a.im);
  }
  for (i = mboProdSpaceDim(h2) / 2; i < mboProdSpaceDim(h2); ++i) {
    EXPECT_DOUBLE_EQ(y[i].re, a.re);
    EXPECT_DOUBLE_EQ(y[i].im, a.im);
  }
  mboTensorOpDestroy(&A);
  mboNumOpDestroy(&Ac);
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
  MboNumOp Cc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(C, &Cc));
  err = mboNumOpMatVec(one, Cc, &x[0], b, &y[0]);
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
    EXPECT_DOUBLE_EQ(y[i].re, expectedResult.re) << "i == " << i;
    EXPECT_DOUBLE_EQ(y[i].im, expectedResult.im) << "i == " << i;
  }
  free(dims);
  mboTensorOpDestroy(&C);
  mboNumOpDestroy(&Cc);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);
}

TEST(MboTensorOp, Flops) {
  MboTensorOp a;
  MboNumOp ac;
  MboElemOp sz;
  MboProdSpace h, h1;
  double flops;

  h = mboProdSpaceCreate(2);
  mboTensorOpNull(h, &a);
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(a, &ac));
  flops = mboNumOpFlops(ac);
  EXPECT_DOUBLE_EQ(flops, 0.0);
  mboTensorOpDestroy(&a);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&ac);

  h = mboProdSpaceCreate(5);
  mboTensorOpIdentity(h, &a);
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(a, &ac));
  flops = mboNumOpFlops(ac);
  EXPECT_DOUBLE_EQ(flops, 5 * 8);
  mboTensorOpDestroy(&a);
  mboNumOpDestroy(&ac);
  mboProdSpaceDestroy(&h);

  h1 = mboProdSpaceCreate(5);
  h = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboProdSpaceMul(h1, &h);
  mboTensorOpIdentity(h, &a);
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(a, &ac));
  flops = mboNumOpFlops(ac);
  EXPECT_DOUBLE_EQ(flops, 5 * 5 * 5 * 8);
  mboTensorOpDestroy(&a);
  mboNumOpDestroy(&ac);
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
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(a, &ac));
  flops = mboNumOpFlops(ac);
  EXPECT_DOUBLE_EQ(flops, (5 * 5 * 5 + 2 * 5 * 5) * 8);
  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&a);
  mboNumOpDestroy(&ac);
  mboProdSpaceDestroy(&h);
  mboProdSpaceDestroy(&h1);
}

TEST(MboTensorOp, MulSelf) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboElemOp eop = mboSigmaZ();
  MboTensorOp A;
  mboTensorOpNull(h, &A);
  mboTensorOpAddTo(eop, 0, A);
  mboTensorOpMul(A, A, &A);

  std::vector<struct MboAmplitude> x(2);
  x[0].re = 134.3;
  x[0].im = 11.0;
  x[1].re = -12.0;
  x[1].im = 0.57;
  std::vector<struct MboAmplitude> y(2);
  struct MboAmplitude one;
  one.re = 1.0;
  one.im = 0.0;
  struct MboAmplitude zero;
  zero.re = 0.0;
  zero.im = 0.0;

  MboNumOp Acomp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(A, &Acomp));
  mboNumOpMatVec(one, Acomp, &x[0], zero, &y[0]);
  EXPECT_DOUBLE_EQ(0.0, y[0].re);
  EXPECT_DOUBLE_EQ(0.0, y[0].im);
  EXPECT_DOUBLE_EQ(2.0 * x[1].re, y[1].re);
  EXPECT_DOUBLE_EQ(2.0 * x[1].im, y[1].im);

  mboNumOpDestroy(&Acomp);
  mboTensorOpDestroy(&A);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, MatVec3BodyOperator) {
  MboGlobInd dim;
  MboProdSpace h1, h2;
  MboTensorOp A, B, C, D;
  struct MboAmplitude a, b, one;
  MboElemOp eop;

  one.re = 1.0;
  one.im = 0.0;
  a = one;
  b.re = 0.0;
  b.im = 0.0;

  // h2 == 2 x 3 x 2 x 2 x 3 x 2
  h1 = mboProdSpaceCreate(2);
  h2 = mboProdSpaceCreate(0);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceDestroy(&h1);
  h1 = mboProdSpaceCreate(3);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceDestroy(&h1);
  h1 = mboProdSpaceCreate(2);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceMul(h2, &h2);
  mboProdSpaceDestroy(&h1);
  dim = mboProdSpaceDim(h2);

  mboTensorOpNull(h2, &D);

  // D == sz1 x sp3 x sm4
  mboProdSpaceDestroy(&h2);
  h2 = mboProdSpaceCreate(3);
  h1 = mboProdSpaceCreate(2);
  mboProdSpaceMul(h1, &h2);
  mboProdSpaceMul(h2, &h1);
  eop =  mboSigmaZ();
  mboTensorOpNull(h1, &A);
  mboTensorOpAddTo(eop, 1, A);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);

  eop = mboSigmaPlus();
  h2 = mboProdSpaceCreate(2);
  mboTensorOpNull(h2, &B);
  mboTensorOpAddTo(eop, 0, B);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h2);

  eop = mboSigmaMinus();
  h1 = mboProdSpaceCreate(3);
  h2 = mboProdSpaceCreate(2);
  mboProdSpaceMul(h1, &h2);
  mboTensorOpNull(h2, &C);
  mboTensorOpAddTo(eop, 0, C);
  mboElemOpDestroy(&eop);
  mboProdSpaceDestroy(&h1);
  mboProdSpaceDestroy(&h2);

  MboTensorOp ops[3] = {A, B, C};
  MBO_STATUS err = mboTensorOpKron(3, ops, &D);
  ASSERT_EQ(MBO_SUCCESS, err);

  MboNumOp Dcomp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(D, &Dcomp));
  std::vector<struct MboAmplitude> x(dim);
  for (int i = 0; i < x.size(); ++i) {
    x[i].re = i;
    x[i].im = -2.7 * i;
  }
  std::vector<struct MboAmplitude> y(dim);
  err = mboNumOpMatVec(a, Dcomp, &x[0], b, &y[0]);
  EXPECT_EQ(MBO_SUCCESS, err);
  EXPECT_DOUBLE_EQ(0.0, y[1].re);
  EXPECT_DOUBLE_EQ(0.0, y[5].re);
  EXPECT_DOUBLE_EQ(-2.0, y[6].re);
  EXPECT_DOUBLE_EQ(-3.0, y[7].re);
  EXPECT_DOUBLE_EQ(0.0, y[8].re);
  EXPECT_DOUBLE_EQ(-14.0, y[18].re);
  EXPECT_DOUBLE_EQ(26.0, y[30].re);

  mboNumOpDestroy(&Dcomp);
  mboTensorOpDestroy(&A);
  mboTensorOpDestroy(&B);
  mboTensorOpDestroy(&C);
  mboTensorOpDestroy(&D);
}

TEST(MboTensorOp, DenseMatrixNull) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboTensorOp null;
  MboNumOp nullc;
  mboTensorOpNull(h, &null);

  MboGlobInd dim = mboProdSpaceDim(h);
  struct MboAmplitude *mat = new struct MboAmplitude[dim * dim];
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(null, &nullc));
  mboNumOpDenseMatrix(nullc, mat);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_DOUBLE_EQ(mat[i].re, 0);
    EXPECT_DOUBLE_EQ(mat[i].im, 0);
  }

  delete [] mat;
  mboTensorOpDestroy(&null);
  mboNumOpDestroy(&nullc);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, DenseMatrixIdentity) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboTensorOp id;
  MboNumOp idc;
  mboTensorOpIdentity(h, &id);

  MboGlobInd dim = mboProdSpaceDim(h);
  struct MboAmplitude *mat = new struct MboAmplitude[dim * dim];
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(id, &idc));
  mboNumOpDenseMatrix(idc, mat);
  for (int i = 0; i < dim * dim; ++i) {
    struct MboAmplitude expectedResult;
    if (i % mboProdSpaceDim(h) == i / mboProdSpaceDim(h)) {
      expectedResult.re = 1.0;
      expectedResult.im = 0.0;
    } else {
      expectedResult.re = 0.0;
      expectedResult.im = 0.0;
    }
    EXPECT_DOUBLE_EQ(mat[i].re, expectedResult.re);
    EXPECT_DOUBLE_EQ(mat[i].im, expectedResult.im);
  }

  delete [] mat;
  mboTensorOpDestroy(&id);
  mboNumOpDestroy(&idc);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsNull) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);
  MboNumOp nullc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(null, &nullc));

  std::vector<MboGlobInd> i(dim + 1);
  mboNumOpRowOffsets(nullc, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(0, i[2]);

  mboTensorOpDestroy(&null);
  mboNumOpDestroy(&nullc);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsIdentity) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);
  MboNumOp idc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(id, &idc));

  std::vector<MboGlobInd> i(dim + 1);
  mboNumOpRowOffsets(idc, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(1, i[1]);
  EXPECT_EQ(2, i[2]);

  mboTensorOpDestroy(&id);
  mboNumOpDestroy(&idc);
  mboProdSpaceDestroy(&h);
}

TEST(MboTensorOp, RowOffsetsIdentityEmptyRange) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);
  std::vector<MboGlobInd> i(dim + 1);
  MboNumOp idc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(id, &idc));
  mboNumOpRowOffsets(idc, 3, 2, &i[0]);
  mboNumOpRowOffsets(idc, 2, 2, &i[0]);
  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&idc);
}

TEST(MboTensorOp, RowOffsetsIdentitySubrange) {
  int dim = 5;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp id;
  mboTensorOpIdentity(h, &id);
  MboNumOp idc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(id, &idc));

  std::vector<MboGlobInd> i(3);
  mboNumOpRowOffsets(idc, 2, 4, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(1, i[1]);
  EXPECT_EQ(2, i[2]);

  mboTensorOpDestroy(&id);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&idc);
}

TEST(MboTensorOp, RowOffsetsSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<MboGlobInd> i(dim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, 2, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(0, i[1]);
  EXPECT_EQ(1, i[2]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
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
  std::vector<MboGlobInd> i(totalDim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, totalDim, &i[0]);
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
  mboNumOpDestroy(&Spc);
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
  std::vector<MboGlobInd> i(totalDim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, totalDim, &i[0]);
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
  mboNumOpDestroy(&Spc);
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
  std::vector<MboGlobInd> i(4);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 5, 8, &i[0]);
  EXPECT_EQ(0, i[0]);
  EXPECT_EQ(2, i[1]);
  EXPECT_EQ(3, i[2]);
  EXPECT_EQ(5, i[3]);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
}

TEST(MboTensorOp, SparseMatrixNull) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<MboGlobInd> i(dim + 1);
  MboNumOp nullc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(null, &nullc));
  mboNumOpRowOffsets(nullc, 0, 2, &i[0]);
  std::vector<MboGlobInd> j(i[dim] + 1);
  std::vector<struct MboAmplitude> a(i[dim] + 1);
  mboNumOpSparseMatrix(nullc, 0, 2, &i[0], &j[0], &a[0]);

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&nullc);
}

TEST(MboTensorOp, SparseMatrixIdentity) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  std::vector<MboGlobInd> i(dim + 1);
  MboNumOp identityc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(identity, &identityc));
  mboNumOpRowOffsets(identityc, 0, 2, &i[0]);
  std::vector<MboGlobInd> j(i[dim]);
  std::vector<struct MboAmplitude> a(i[dim]);
  mboNumOpSparseMatrix(identityc, 0, 2, &i[0], &j[0], &a[0]);
  EXPECT_EQ(0, j[0]);
  EXPECT_EQ(1, j[1]);
  EXPECT_DOUBLE_EQ(1.0, a[0].re);
  EXPECT_DOUBLE_EQ(0.0, a[0].im);
  EXPECT_DOUBLE_EQ(1.0, a[1].re);
  EXPECT_DOUBLE_EQ(0.0, a[1].im);

  mboTensorOpDestroy(&identity);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&identityc);
}

TEST(MboTensorOp, SparseMatrixSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<MboGlobInd> i(dim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, 2, &i[0]);
  std::vector<MboGlobInd> j(i[dim]);
  std::vector<struct MboAmplitude> a(i[dim]);
  mboNumOpSparseMatrix(Spc, 0, 2, &i[0], &j[0], &a[0]);
  EXPECT_EQ(0, j[0]);
  EXPECT_DOUBLE_EQ(1.0, a[0].re);
  EXPECT_DOUBLE_EQ(0.0, a[0].im);

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
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
  std::vector<MboGlobInd> i(totalDim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, totalDim, &i[0]);
  std::vector<MboGlobInd> j(i[mboProdSpaceDim(h)]);
  std::vector<struct MboAmplitude> a(i[mboProdSpaceDim(h)]);
  mboNumOpSparseMatrix(Spc, 0, mboProdSpaceDim(h), &i[0], &j[0], &a[0]);
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
	  EXPECT_DOUBLE_EQ(1.0, a[i].re) << "i = " << i;
	  EXPECT_DOUBLE_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
}

static void printMatrixPattern(MboNumOp op) {
  MboGlobInd dim = mboProdSpaceDim(mboNumOpGetSpace(op));
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboNumOpDenseMatrix(op, &mat[0]);
  std::cout << "\n\n";
  assert(dim > 0);
  for (size_t i = 0; i < (size_t)dim; ++i) {
    for (size_t j = 0; j < (size_t)dim; ++j) {
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
  std::vector<MboGlobInd> i(totalDim + 1);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, 0, totalDim, &i[0]);
  std::vector<MboGlobInd> j(i[mboProdSpaceDim(h)]);
  std::vector<struct MboAmplitude> a(i[mboProdSpaceDim(h)]);
  mboNumOpSparseMatrix(Spc, 0, mboProdSpaceDim(h), &i[0], &j[0], &a[0]);
  for (int m = 0; m < totalDim; ++m) {
    std::sort(&j[0] + i[m], &j[0] + i[m + 1]);
  }
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
    EXPECT_DOUBLE_EQ(1.0, a[i].re) << "i = " << i;
    EXPECT_DOUBLE_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
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
  std::vector<MboGlobInd> i(totalDim + 1);
  MboGlobInd rmin = 3;
  MboGlobInd rmax = 11;
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpRowOffsets(Spc, rmin, rmax, &i[0]);
  std::vector<MboGlobInd> j(i[rmax - rmin]);
  std::vector<struct MboAmplitude> a(i[rmax - rmin]);
  mboNumOpSparseMatrix(Spc, rmin, rmax, &i[0], &j[0], &a[0]);
  for (int m = 0; m < rmax - rmin; ++m) {
    std::sort(&j[0] + i[m], &j[0] + i[m + 1]);
  }
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
    EXPECT_DOUBLE_EQ(1.0, a[i].re) << "i = " << i;
    EXPECT_DOUBLE_EQ(0.0, a[i].im) << "i = " << i;
  }

  mboElemOpDestroy(&sp);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Spc);
}

TEST(MboNumOp, GetDiagonalNull) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  std::vector<struct MboAmplitude> diag(dim);
  MboNumOp nullc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(null, &nullc));
  mboNumOpDiagonal(nullc, 0, dim, &diag[0]);

  mboTensorOpDestroy(&null);
  mboNumOpDestroy(&nullc);
  mboProdSpaceDestroy(&h);
}

TEST(MboNumOp, GetDiagonalIdentity) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  std::vector<struct MboAmplitude> diag(dim);
  MboNumOp identityc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(identity, &identityc));
  mboNumOpDiagonal(identityc, 0, dim, &diag[0]);
  for (int i = 0; i < dim; ++i) {
    EXPECT_DOUBLE_EQ(1.0, diag[i].re) << "i = " << i;
    EXPECT_DOUBLE_EQ(0.0, diag[i].im) << "i = " << i;
  }

  mboTensorOpDestroy(&identity);
  mboNumOpDestroy(&identityc);
  mboProdSpaceDestroy(&h);
}

TEST(MboNumOp, GetDiagonalSigmaPlus) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 0, Sp);

  std::vector<struct MboAmplitude> diag(dim);
  MboNumOp Spc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sp, &Spc));
  mboNumOpDiagonal(Spc, 0, dim, &diag[0]);
  EXPECT_DOUBLE_EQ(0.0, diag[0].re);
  EXPECT_DOUBLE_EQ(0.0, diag[0].im);
  EXPECT_DOUBLE_EQ(0.0, diag[1].re);
  EXPECT_DOUBLE_EQ(0.0, diag[1].im);

  mboElemOpDestroy(&sp);
  mboNumOpDestroy(&Spc);
  mboTensorOpDestroy(&Sp);
  mboProdSpaceDestroy(&h);
}

TEST(MboNumOp, GetDiagonalSigmaZ) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp Sz;
  mboTensorOpNull(h, &Sz);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 0, Sz);

  std::vector<struct MboAmplitude> diag(dim);
  MboNumOp Szc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sz, &Szc));
  mboNumOpDiagonal(Szc, 0, dim, &diag[0]);
  EXPECT_DOUBLE_EQ(-1.0, diag[0].re);
  EXPECT_DOUBLE_EQ(0.0, diag[0].im);
  EXPECT_DOUBLE_EQ(1.0, diag[1].re);
  EXPECT_DOUBLE_EQ(0.0, diag[1].im);

  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&Sz);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Szc);
}

TEST(MboNumOp, GetDiagonalSigmaZInBiggerSpace) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Sz;
  mboTensorOpNull(h, &Sz);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 2, Sz);

  std::vector<struct MboAmplitude> diag(mboProdSpaceDim(h));
  MboNumOp Szc;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Sz, &Szc));
  mboNumOpDiagonal(Szc, 0, mboProdSpaceDim(h), &diag[0]);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (j / 2 == 0) {
        EXPECT_DOUBLE_EQ(-1.0, diag[i * 4 + j].re) << " i == " << i
                                                  << " j == " << j;
        EXPECT_DOUBLE_EQ(0.0, diag[i * 4 + j].im) << " i == " << i
                                                 << " j == " << j;
      } else {
        EXPECT_DOUBLE_EQ(1.0, diag[i * 4 + j].re) << " i == " << i
                                                 << " j == " << j;
        EXPECT_DOUBLE_EQ(0.0, diag[i * 4 + j].im) << " i == " << i
                                                 << " j == " << j;
      }
    }
  }

  mboElemOpDestroy(&sz);
  mboNumOpDestroy(&Szc);
  mboTensorOpDestroy(&Sz);
  mboProdSpaceDestroy(&h);
}

TEST(MboNumOp, DeleteDiagonalNull) {
  int dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp null;
  mboTensorOpNull(h, &null);

  MboNumOp null_comp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(null, &null_comp));
  mboNumOpDeleteDiagonal(null_comp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboNumOpDenseMatrix(null_comp, &mat[0]);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_DOUBLE_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(0, mat[i].im) << " i == " << i;
  }

  mboTensorOpDestroy(&null);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&null_comp);
}

TEST(MboNumOp, DeleteDiagonalIdentity) {
  int dim = 3;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboTensorOp identity;
  mboTensorOpIdentity(h, &identity);

  MboNumOp identity_comp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(identity, &identity_comp));
  mboNumOpDeleteDiagonal(identity_comp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboNumOpDenseMatrix(identity_comp, &mat[0]);
  for (int i = 0; i < dim * dim; ++i) {
    EXPECT_DOUBLE_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(0, mat[i].im) << " i == " << i;
  }

  mboTensorOpDestroy(&identity);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&identity_comp);
}

TEST(MboNumOp, DeleteDiagonalSz) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Op;
  mboTensorOpNull(h, &Op);
  MboElemOp sz = mboSigmaZ();
  mboTensorOpAddTo(sz, 1, Op);

  MboNumOp Op_comp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Op, &Op_comp));
  mboNumOpDeleteDiagonal(Op_comp);
  MboGlobInd D = mboProdSpaceDim(h);
  std::vector<struct MboAmplitude> mat(D * D);
  mboNumOpDenseMatrix(Op_comp, &mat[0]);
  for (int i = 0; i < D * D; ++i) {
    EXPECT_DOUBLE_EQ(0, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(0, mat[i].im) << " i == " << i;
  }

  mboElemOpDestroy(&sz);
  mboTensorOpDestroy(&Op);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Op_comp);
}

TEST(MboNumOp, DeleteDiagonalSp) {
  int dim = 2;
  MboProdSpace h = mboProdSpaceCreate(dim);
  mboProdSpaceMul(h, &h);
  mboProdSpaceMul(h, &h);
  MboTensorOp Op;
  mboTensorOpNull(h, &Op);
  MboElemOp sp = mboSigmaPlus();
  mboTensorOpAddTo(sp, 1, Op);

  MboNumOp Op_comp;
  ASSERT_EQ(MBO_SUCCESS, mboNumOpCompile(Op, &Op_comp));
  mboNumOpDeleteDiagonal(Op_comp);
  MboGlobInd D = mboProdSpaceDim(h);
  std::vector<MboGlobInd> I(D + 1);
  mboNumOpRowOffsets(Op_comp, 0, D, &I[0]);
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
  mboNumOpDestroy(&Op_comp);
}
