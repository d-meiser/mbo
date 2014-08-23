#include <gtest/gtest.h>
#include <MboVec.h>
#include <MboAmplitude.h>

TEST(MboVec, Create) {
  MboVec v;
  mboVecCreate(5, &v);
  EXPECT_NE(v, (MboVec)0);
  mboVecDestroy(&v);
}

TEST(MboVec, CreateZeroDimension) {
  MboVec v = 0;
  MBO_STATUS err = mboVecCreate(0, &v);
  EXPECT_EQ(err, MBO_ILLEGAL_DIMENSION);
  mboVecDestroy(&v);
}

TEST(MboVec, CreateNegativeDimension) {
  MboVec v = 0;
  MBO_STATUS err = mboVecCreate(-10, &v);
  EXPECT_EQ(err, MBO_ILLEGAL_DIMENSION);
  mboVecDestroy(&v);
}

TEST(MboVec, Destroy) {
  MboVec v;
  mboVecCreate(5, &v);
  mboVecDestroy(&v);
  EXPECT_EQ(v, (MboVec)0);
}

TEST(MboVec, CheckPass) {
  MboVec v;
  mboVecCreate(5, &v);
  EXPECT_EQ(mboVecCheck(v), 0);
  mboVecDestroy(&v);
}

TEST(MboVec, GetViewRW) {
  MboVec v;
  struct MboAmplitude *arr;
  int n = 5;
  mboVecCreate(n, &v);
  int err = mboVecGetViewRW(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_NE(arr, (struct MboAmplitude *)0);
  for (int i = 0; i < n; ++i) {
    arr[i].re = i;
    arr[i].im = -2.0 * i;
  }
  err = mboVecReleaseView(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_EQ(arr, (struct MboAmplitude *)0);
  mboVecDestroy(&v);
}

TEST(MboVec, GetViewR) {
  int n = 5;
  MboVec v;
  struct MboAmplitude *arr;
  mboVecCreate(n, &v);
  int err = mboVecGetViewR(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_NE(arr, (struct MboAmplitude *)0);
  err = mboVecGetViewR(v, &arr);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  EXPECT_EQ(arr, (struct MboAmplitude *)0);
  err = mboVecReleaseView(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);

  err = mboVecGetViewRW(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (int i = 0; i < n; ++i) {
    arr[i].re = i;
    arr[i].im = -2.0 * i;
  }
  mboVecReleaseView(v, &arr);
  err = mboVecGetViewR(v, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (int i = 0; i < n; ++i) {
    EXPECT_FLOAT_EQ(arr[i].re, i);
    EXPECT_FLOAT_EQ(arr[i].im, -2.0 * i);
  }
  mboVecReleaseView(v, &arr);
  mboVecDestroy(&v);
}

TEST(MboVec, AXPY) {
  MboVec x, y;
  struct MboAmplitude a;
  a.re = 3.7;
  a.im = 2.0;
  int d = 2;
  mboVecCreate(d, &x);
  mboVecCreate(d, &y);
  MBO_STATUS err = mboVecAXPY(&a, x, y);
  EXPECT_EQ(err, MBO_SUCCESS);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
}

TEST(MboVec, AXPYInUse) {
  MboVec x, y;
  struct MboAmplitude a, *arr;
  a.re = 3.7;
  a.im = 2.0;
  int d = 2;
  mboVecCreate(d, &x);
  mboVecCreate(d, &y);
  mboVecGetViewR(x, &arr);
  MBO_STATUS err = mboVecAXPY(&a, x, y);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
}

TEST(MboVec, AXPYDimensionsWrong) {
  MboVec x, y;
  struct MboAmplitude a;
  MBO_STATUS err;
  a.re = 3.7;
  a.im = 2.0;
  int d = 2;
  mboVecCreate(d + 1, &x);
  mboVecCreate(d, &y);
  err = mboVecAXPY(&a, x, y);
  EXPECT_EQ(err, MBO_DIMENSIONS_MISMATCH);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
}

TEST(MboVec, Set) {
  int d = 2;
  MboVec x;
  struct MboAmplitude a, *arr;
  a.re = 3.7;
  a.im = 2.0;
  mboVecCreate(d, &x);
  MBO_STATUS err = mboVecSet(&a, x);
  EXPECT_EQ(err, MBO_SUCCESS);
  mboVecGetViewR(x, &arr);
  err = mboVecSet(&a, x);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  err = mboVecReleaseView(x, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecDestroy(&x);
  EXPECT_EQ(err, MBO_SUCCESS);
}

TEST(MboVec, SetRandom) {
  int d = 2;
  MboVec x;
  struct MboAmplitude *arr, zero;
  MBO_STATUS err = mboVecCreate(d, &x);
  ASSERT_EQ(err, MBO_SUCCESS);
  mboVecGetViewR(x, &arr);
  err = mboVecSetRandom(x);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  err = mboVecReleaseView(x, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecDestroy(&x);
  EXPECT_EQ(err, MBO_SUCCESS);

  err = mboVecCreate(d, &x);
  ASSERT_EQ(err, MBO_SUCCESS);
  zero.re = 0;
  zero.im = 0;
  err = mboVecSet(&zero, x);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecSetRandom(x);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecGetViewR(x, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  for (int i = 0; i < d; ++i) {
    EXPECT_NE(arr[i].re, 0);
    EXPECT_NE(arr[i].im, 0);
  }
  err = mboVecReleaseView(x, &arr);
  EXPECT_EQ(err, MBO_SUCCESS);
  err = mboVecDestroy(&x);
  EXPECT_EQ(err, MBO_SUCCESS);
}

TEST(MboVec, Duplicate) {
  MboVec x, y;
  int d = 5;
  mboVecCreate(d, &x);
  mboVecSetRandom(x);
  MBO_STATUS err = mboVecDuplicate(x, &y);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_EQ(mboVecDim(x), mboVecDim(y));
  struct MboAmplitude *xarr, *yarr;
  mboVecGetViewR(x, &xarr);
  mboVecGetViewR(y, &yarr);
  for (int i = 0; i < mboVecDim(x); ++i) {
    EXPECT_FLOAT_EQ(xarr[i].re, yarr[i].re);
    EXPECT_FLOAT_EQ(xarr[i].im, yarr[i].im);
  }
  mboVecDestroy(&y);
  mboVecDestroy(&x);
}

TEST(MboVec, DuplicateInUse) {
  MboVec x = 0, y = 0;
  struct MboAmplitude *arr;
  int d = 5;
  mboVecCreate(d, &x);
  mboVecSetRandom(x);
  mboVecGetViewR(x, &arr);
  MBO_STATUS err = mboVecDuplicate(x, &y);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  mboVecDestroy(&y);
  mboVecDestroy(&x);
}

static void buildArrays(int n, MboLocInd *dims, struct MboAmplitude ***arrays) {
  int i, j;
  *arrays = (struct MboAmplitude **)malloc(n * sizeof(**arrays));
  for (i = 0; i < n; ++i) {
    (*arrays)[i] = (struct MboAmplitude *)malloc(dims[i] * sizeof(***arrays));
    for (j = 0; j < dims[i]; ++j) {
      (*arrays)[i][j].re = 2.0 * i - 3.0 * j;
      (*arrays)[i][j].im = -2.1 * i + 3.1 * j;
    }
  }
}

static void freeArrays(int n, struct MboAmplitude ***arrays) {
  int i;
  for (i = 0; i < n; ++i) {
    free((*arrays)[i]);
  }
  free(*arrays);
  *arrays = 0;
}

TEST(MboVec, KronDimensionsMismatch) {
  struct MboAmplitude **arrays;
  struct MboAmplitude zero, tmp;
  MboVec x;
  zero.re = 0;
  zero.im = 0;
  int n = 1;
  MboLocInd dims[1];
  dims[0] = 2;
  buildArrays(n, dims, &arrays);
  mboVecCreate(dims[0] + 1, &x);
  MBO_STATUS err = mboVecKron(n, dims, arrays, x);
  EXPECT_EQ(err, MBO_DIMENSIONS_MISMATCH);
  mboVecDestroy(&x);
  freeArrays(n, &arrays);
}

TEST(MboVec, KronOne) {
  int n = 1;
  MboLocInd dims[1];
  dims[0] = 2;
  struct MboAmplitude **arrays;
  buildArrays(n, dims, &arrays);
  MboVec x;
  mboVecCreate(dims[0], &x);
  struct MboAmplitude zero;
  zero.re = 0;
  zero.im = 0;
  mboVecSet(&zero, x);
  MBO_STATUS err = mboVecKron(n, dims, arrays, x);
  ASSERT_EQ(err, MBO_SUCCESS);
  struct MboAmplitude *arr;
  mboVecGetViewR(x, &arr);
  EXPECT_FLOAT_EQ(arr[0].re, arrays[0][0].re);
  EXPECT_FLOAT_EQ(arr[0].im, arrays[0][0].im);
  EXPECT_FLOAT_EQ(arr[1].re, arrays[0][1].re);
  EXPECT_FLOAT_EQ(arr[1].im, arrays[0][1].im);
  mboVecReleaseView(x, &arr);
  mboVecDestroy(&x);
  freeArrays(n, &arrays);
}

TEST(MboVec, KronTwo) {
  int n = 2;
  int dims[2];
  dims[0] = 2;
  dims[1] = 3;
  struct MboAmplitude **arrays;
  buildArrays(n, dims, &arrays);
  MboVec x;
  mboVecCreate(dims[0] * dims[1], &x);
  struct MboAmplitude zero;
  zero.re = 0;
  zero.im = 0;
  mboVecSet(&zero, x);
  mboVecKron(n, dims, arrays, x);
  struct MboAmplitude *arr;
  mboVecGetViewR(x, &arr);
  for (int i = 0; i < dims[0]; ++i) {
    for (int j = 0; j < dims[1]; ++j) {
      struct MboAmplitude tmp;
      tmp.re =
          arrays[0][i].re * arrays[1][j].re - arrays[0][i].im * arrays[1][j].im;
      tmp.im =
          arrays[0][i].re * arrays[1][j].im + arrays[0][i].im * arrays[1][j].re;
      EXPECT_FLOAT_EQ(tmp.re, arr[i * dims[1] + j].re);
      EXPECT_FLOAT_EQ(tmp.im, arr[i * dims[1] + j].im);
    }
  }
  mboVecReleaseView(x, &arr);
  mboVecDestroy(&x);
  freeArrays(n, &arrays);
}

TEST(MboVec, Dim) {
  MboGlobInd d;
  MboVec x;
  mboVecCreate(10l, &x);
  d = mboVecDim(x);
  EXPECT_EQ(d, 10);
  mboVecDestroy(&x);
}

static void f1(int n, MboLocInd *dims, MboLocInd *indices, void *ctx,
               struct MboAmplitude *x) {
  int i;
  double value = 0.0, b = 1.0;
  for (i = n - 1; i >= 0; --i) {
    value += indices[i] * b;
    b *= 10.0;
  }
  x->re = value;
  x->im = 0;
}

static void f2(int n, MboLocInd *dims, MboLocInd *indices, void *ctx,
               struct MboAmplitude *x) {
  int *i = (int *)ctx;
  x->re = *i;
  x->im = 0;
  ++*i;
}

TEST(MboVec, Map) {
  MboLocInd dims[] = {2, 3, 2};
  MboVec x;
  struct MboAmplitude *xarr;

  mboVecCreate(12, &x);
  mboVecMap(3, dims, f1, 0, x);
  mboVecGetViewR(x, &xarr);
  EXPECT_FLOAT_EQ(xarr[0].re, 0);
  EXPECT_FLOAT_EQ(xarr[0].im, 0);
  EXPECT_FLOAT_EQ(xarr[1].re, 1);
  EXPECT_FLOAT_EQ(xarr[2].re, 10);
  EXPECT_FLOAT_EQ(xarr[3].re, 11);
  EXPECT_FLOAT_EQ(xarr[4].re, 20);
  EXPECT_FLOAT_EQ(xarr[5].re, 21);
  EXPECT_FLOAT_EQ(xarr[6].re, 100);
  EXPECT_FLOAT_EQ(xarr[7].re, 101);
  EXPECT_FLOAT_EQ(xarr[8].re, 110);
  EXPECT_FLOAT_EQ(xarr[9].re, 111);
  EXPECT_FLOAT_EQ(xarr[10].re, 120);
  EXPECT_FLOAT_EQ(xarr[11].re, 121);
  mboVecReleaseView(x, &xarr);
  mboVecDestroy(&x);

  mboVecCreate(12, &x);
  int i = 0;
  mboVecMap(3, dims, f2, &i, x);
  mboVecGetViewR(x, &xarr);
  EXPECT_EQ(i, 12);
  for (i = 0; i < 12; ++i) {
    EXPECT_FLOAT_EQ(xarr[i].re, i);
    EXPECT_FLOAT_EQ(xarr[i].im, 0);
  }
  mboVecReleaseView(x, &xarr);
  mboVecDestroy(&x);
}

TEST(MboVec, Dot) {
  MboVec x, y;
  struct MboAmplitude result, *dummy;
  mboVecCreate(3, &x);
  mboVecCreate(2, &y);
  MBO_STATUS err = mboVecDot(x, y, &result);
  EXPECT_EQ(err, MBO_DIMENSIONS_MISMATCH);
  mboVecDestroy(&x);
  mboVecDestroy(&y);

  mboVecCreate(2, &x);
  mboVecCreate(2, &y);
  mboVecGetViewRW(x, &dummy);
  err = mboVecDot(x, y, &result);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  mboVecDestroy(&x);
  mboVecDestroy(&y);

  mboVecCreate(2, &x);
  mboVecGetViewRW(x, &dummy);
  dummy[0].re = 2.0;
  dummy[0].im = 1.0;
  dummy[1].re = 3.0;
  dummy[1].im = -1.0;
  mboVecReleaseView(x, &dummy);
  mboVecCreate(2, &y);
  mboVecGetViewRW(y, &dummy);
  dummy[0].re = 1.0;
  dummy[0].im = 3.0;
  dummy[1].re = 10.0;
  dummy[1].im = 15.0;
  mboVecReleaseView(y, &dummy);
  err = mboVecDot(x, y, &result);
  EXPECT_EQ(err, MBO_SUCCESS);
  EXPECT_FLOAT_EQ(result.re,
                  2.0 * 1.0 + 1.0 * 3.0 + 3.0 * 10.0 + (-1.0) * 15.0);
  EXPECT_FLOAT_EQ(result.im,
                  2.0 * 3.0 - 1.0 * 1.0 + 3.0 * 15.0 - (-1.0) * 10.0);
  mboVecDestroy(&x);
  mboVecDestroy(&y);
}

TEST(MboVec, UnitVector) {
  MboVec x;
  mboVecCreate(2, &x);
  MBO_STATUS err = mboVecUnitVector(2, x);
  EXPECT_EQ(err, MBO_ILLEGAL_DIMENSION);
  mboVecDestroy(&x);

  mboVecCreate(2, &x);
  struct MboAmplitude *xarr;
  mboVecGetViewRW(x, &xarr);
  err = mboVecUnitVector(1, x);
  EXPECT_EQ(err, MBO_VEC_IN_USE);
  mboVecReleaseView(x, &xarr);
  mboVecDestroy(&x);

  mboVecCreate(2, &x);
  err = mboVecUnitVector(1, x);
  EXPECT_EQ(err, MBO_SUCCESS);
  mboVecGetViewRW(x, &xarr);
  EXPECT_FLOAT_EQ(xarr[0].re, 0);
  EXPECT_FLOAT_EQ(xarr[0].im, 0);
  EXPECT_FLOAT_EQ(xarr[1].re, 1);
  EXPECT_FLOAT_EQ(xarr[1].im, 0);
  mboVecReleaseView(x, &xarr);
  mboVecDestroy(&x);
}
