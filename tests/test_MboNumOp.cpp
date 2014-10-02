#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
#include <MboAmplitude.h>
#include <vector>


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

