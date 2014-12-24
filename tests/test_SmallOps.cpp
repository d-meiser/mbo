#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
#include <MboAmplitude.h>
#include <MboElemOp.h>
#include <cmath>

TEST(SmallOps, SigmaP) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboElemOp sp = mboSigmaPlus();
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(4);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(4);
  expectedResult[0].re = 0;
  expectedResult[0].im = 0;
  expectedResult[1].re = 0;
  expectedResult[1].im = 0;
  expectedResult[2].re = 1;
  expectedResult[2].im = 0;
  expectedResult[3].re = 0;
  expectedResult[3].im = 0;
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, SigmaM) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboElemOp sp = mboSigmaMinus();
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(4);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(4);
  expectedResult[0].re = 0;
  expectedResult[0].im = 0;
  expectedResult[1].re = 1;
  expectedResult[1].im = 0;
  expectedResult[2].re = 0;
  expectedResult[2].im = 0;
  expectedResult[3].re = 0;
  expectedResult[3].im = 0;
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, SigmaZ) {
  MboProdSpace h = mboProdSpaceCreate(2);
  MboElemOp sp = mboSigmaZ();
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(4);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(4);
  expectedResult[0].re = -1;
  expectedResult[0].im = 0;
  expectedResult[1].re = 0;
  expectedResult[1].im = 0;
  expectedResult[2].re = 0;
  expectedResult[2].im = 0;
  expectedResult[3].re = 1;
  expectedResult[3].im = 0;
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, Eye) {
  MboLocInd dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboElemOp sp = mboEye(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(dim * dim);
  for (int i = 0; i < dim * dim; ++i) {
    if ((i % dim) == (i / dim)) {
      expectedResult[i].re = 1.0;
    } else {
      expectedResult[i].re = 0.0;
    }
    expectedResult[i].im = 0.0;
  }
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, NumOp) {
  MboLocInd dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboElemOp sp = mboNumOp(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(dim * dim);
  for (int i = 0; i < dim * dim; ++i) {
    int row = i / dim;
    int col = i % dim;
    if (col == row) {
      expectedResult[i].re = row;
    } else {
      expectedResult[i].re = 0.0;
    }
    expectedResult[i].im = 0.0;
  }
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, AnnihilationOp) {
  MboLocInd dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboElemOp sp = mboAnnihilationOp(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(dim * dim);
  for (int i = 0; i < dim * dim; ++i) {
    int row = i / dim;
    int col = i % dim;
    if (col == row + 1) {
      expectedResult[i].re = sqrt(col);
    } else {
      expectedResult[i].re = 0.0;
    }
    expectedResult[i].im = 0.0;
  }
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}

TEST(SmallOps, CreationOp) {
  MboLocInd dim = 4;
  MboProdSpace h = mboProdSpaceCreate(dim);
  MboElemOp sp = mboCreationOp(dim);
  MboTensorOp Sp;
  mboTensorOpNull(h, &Sp);
  mboTensorOpAddTo(sp, 0, Sp);
  std::vector<struct MboAmplitude> mat(dim * dim);
  MboNumOp Sp_comp = mboNumOpCompile(Sp);
  mboNumOpDenseMatrix(Sp_comp, &mat[0]);
  std::vector<struct MboAmplitude> expectedResult(dim * dim);
  for (int i = 0; i < dim * dim; ++i) {
    int row = i / dim;
    int col = i % dim;
    if (col + 1 == row) {
      expectedResult[i].re = sqrt(row);
    } else {
      expectedResult[i].re = 0.0;
    }
    expectedResult[i].im = 0.0;
  }
  for (size_t i = 0; i < mat.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedResult[i].re, mat[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult[i].im, mat[i].im) << " i == " << i;
  }
  mboTensorOpDestroy(&Sp);
  mboElemOpDestroy(&sp);
  mboProdSpaceDestroy(&h);
  mboNumOpDestroy(&Sp_comp);
}
