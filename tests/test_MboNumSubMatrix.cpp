#include <gtest/gtest.h>
#include <MboNumSubMatrix.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
#include <MboProdSpace.h>
#include <MboAmplitude.h>

class MboNumSubMatrixFixture : public ::testing::Test {
 public:
  MboNumOp Jx;
  MboProdSpace h;
  static const int n = 4;

  void SetUp() {
    h = buildSpace(n);
    Jx = buildJx(h);
  }

  void TearDown() {
    mboProdSpaceDestroy(&h);
    mboNumOpDestroy(&Jx);
  }

 private:
  MboProdSpace buildSpace(int n) {
    MboProdSpace h1, hTot;
    int i;

    h1 = mboProdSpaceCreate(2);
    hTot = mboProdSpaceCreate(0);
    for (i = 0; i < n; ++i) {
      mboProdSpaceMul(h1, &hTot);
    }
    mboProdSpaceDestroy(&h1);
    return hTot;
  }
  MboNumOp buildJx(MboProdSpace h) {
    MboElemOp sp, sm;
    MboTensorOp Jx;
    MboNumOp Jx_compiled;
    struct MboAmplitude pointFive;
    int i;

    sp = mboSigmaPlus();
    sm = mboSigmaMinus();
    pointFive.re = 0.5;
    pointFive.im = 0.0;
    mboTensorOpNull(h, &Jx);
    for (i = 0; i < mboProdSpaceSize(h); ++i) {
      mboTensorOpAddScaledTo(&pointFive, sm, i, Jx);
      mboTensorOpAddScaledTo(&pointFive, sp, i, Jx);
    }
    mboElemOpDestroy(&sp);
    mboElemOpDestroy(&sm);
    Jx_compiled = mboNumOpCompile(Jx);
    mboTensorOpDestroy(&Jx);
    return Jx_compiled;
  }
};

TEST_F(MboNumSubMatrixFixture, Create) {
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 10, 12, 13, 17);
  ASSERT_NE((void*)0, m);
  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, SetTile) {
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 10, 12, 13, 17);
  mboNumSubMatrixSetTile(m, 12, 13, 20, 30);
  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, MatVec) {
  MboGlobInd dim = mboProdSpaceDim(h);
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 0, dim, 0, dim);

  std::vector<struct MboAmplitude> x(dim);
  struct MboAmplitude one = {1, 0};
  std::fill(x.begin(), x.end(), one);
  std::vector<struct MboAmplitude> y(dim);
  std::fill(y.begin(), y.end(), one);

  struct MboAmplitude zero = {0, 0};
  mboNumSubMatrixMatVec(one, m, &x[0], zero, &y[0]);
  for (MboGlobInd i = 0; i < dim; ++i) {
	  EXPECT_FLOAT_EQ(2.0, y[i].re) << "i == " << i;
	  EXPECT_FLOAT_EQ(0.0, y[i].im) << "i == " << i;
  }

  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, MatVecFirstRow) {
  MboGlobInd dim = mboProdSpaceDim(h);
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 0, 1, 0, dim);

  std::vector<struct MboAmplitude> x(dim);
  struct MboAmplitude one = {1, 0};
  std::fill(x.begin(), x.end(), one);
  std::vector<struct MboAmplitude> y(dim);
  std::fill(y.begin(), y.end(), one);

  struct MboAmplitude zero = {0, 0};
  mboNumSubMatrixMatVec(one, m, &x[0], zero, &y[0]);
  EXPECT_FLOAT_EQ(2.0, y[0].re) << "i == " << 0;
  EXPECT_FLOAT_EQ(0.0, y[0].im) << "i == " << 0;
  for (MboGlobInd i = 1; i < dim; ++i) {
	  EXPECT_FLOAT_EQ(1.0, y[i].re) << "i == " << i;
	  EXPECT_FLOAT_EQ(0.0, y[i].im) << "i == " << i;
  }

  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, MatVecFirstColumn) {
  MboGlobInd dim = mboProdSpaceDim(h);
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 0, dim, 0, 1);

  std::vector<struct MboAmplitude> x(dim);
  struct MboAmplitude one = {1, 0};
  std::fill(x.begin(), x.end(), one);
  std::vector<struct MboAmplitude> y(dim);
  std::fill(y.begin(), y.end(), one);

  struct MboAmplitude zero = {0, 0};
  mboNumSubMatrixMatVec(one, m, &x[0], zero, &y[0]);
  double expectedRealParts[] = {0,   0.5, 0.5, 0, 0.5, 0, 0, 0,
                                0.5, 0,   0,   0, 0,   0, 0, 0};
  for (MboGlobInd i = 1; i < dim; ++i) {
    EXPECT_FLOAT_EQ(expectedRealParts[i], y[i].re) << "i == " << i;
    EXPECT_FLOAT_EQ(0, y[i].im) << "i == " << i;
  }

  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, MatVecTopLeftCorner) {
  MboGlobInd dim = mboProdSpaceDim(h);
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 0, 2, 0, 2);

  std::vector<struct MboAmplitude> x(2);
  struct MboAmplitude one = {1, 0};
  std::fill(x.begin(), x.end(), one);
  std::vector<struct MboAmplitude> y(2);
  std::fill(y.begin(), y.end(), one);

  struct MboAmplitude zero = {0, 0};
  mboNumSubMatrixMatVec(one, m, &x[0], zero, &y[0]);
  EXPECT_FLOAT_EQ(0.5, y[0].re);
  EXPECT_FLOAT_EQ(0.0, y[0].im);
  EXPECT_FLOAT_EQ(0.5, y[1].re);
  EXPECT_FLOAT_EQ(0.0, y[1].im);

  mboNumSubMatrixDestroy(&m);
}

TEST_F(MboNumSubMatrixFixture, MatVecTopRightCorner) {
  MboGlobInd dim = mboProdSpaceDim(h);
  std::cout << "dim == " << dim << std::endl;
  MboNumSubMatrix m = mboNumSubMatrixCreate(Jx, 0, 2, 14, 16);

  std::vector<struct MboAmplitude> x(2);
  struct MboAmplitude one = {1, 0};
  std::fill(x.begin(), x.end(), one);
  std::vector<struct MboAmplitude> y(2);
  std::fill(y.begin(), y.end(), one);

  struct MboAmplitude zero = {0, 0};
  mboNumSubMatrixMatVec(one, m, &x[0], zero, &y[0]);
  EXPECT_FLOAT_EQ(0.0, y[0].re);
  EXPECT_FLOAT_EQ(0.0, y[0].im);
  EXPECT_FLOAT_EQ(0.0, y[1].re);
  EXPECT_FLOAT_EQ(0.0, y[1].im);

  mboNumSubMatrixDestroy(&m);
}
