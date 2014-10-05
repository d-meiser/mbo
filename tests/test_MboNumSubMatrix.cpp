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
  ASSERT_NE((void*)0, m);
  mboNumSubMatrixSetTile(m, 12, 13, 20, 30);
  mboNumSubMatrixDestroy(&m);
}

