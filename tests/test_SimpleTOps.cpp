#include <gtest/gtest.h>
#include <SimpleTOp.h>
#include <vector>

TEST(SimpleTOp, Kron) {
  struct SimpleTOp a, b, c;

  a.numFactors = 0;
  a.embeddings = 0;
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 0, &b, &c);
  EXPECT_EQ(c.numFactors, 0);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 0;
  a.embeddings = 0;
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 1;
  c.embeddings = (struct Embedding*)malloc(sizeof(*c.embeddings));
  c.embeddings[0].i = 0;
  mboElemOpCreate(&c.embeddings[0].op);
  kronSimpleTOps(&a, 0, &b, &c);
  EXPECT_EQ(c.numFactors, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 1;
  a.embeddings = (struct Embedding*)malloc(sizeof(*a.embeddings));
  a.embeddings[0].i = 0;
  mboElemOpCreate(&a.embeddings[0].op);
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 1, &b, &c);
  EXPECT_EQ(c.numFactors, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 1;
  a.embeddings = (struct Embedding*)malloc(sizeof(*a.embeddings));
  a.embeddings[0].i = 0;
  mboElemOpCreate(&a.embeddings[0].op);
  b.numFactors = 1;
  b.embeddings = (struct Embedding*)malloc(sizeof(*b.embeddings));
  b.embeddings[0].i = 0;
  mboElemOpCreate(&b.embeddings[0].op);
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 1, &b, &c);
  EXPECT_EQ(c.numFactors, 2);
  EXPECT_EQ(c.embeddings[0].i, 0);
  EXPECT_EQ(c.embeddings[1].i, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);
}

TEST(SimpleTOp, DenseMatrixIdentityOneSpace) {
  MboLocInd d = 2;
  MboProdSpace h = mboProdSpaceCreate(d);
  struct SimpleTOp sto;
  sto.numFactors = 0;
  sto.embeddings = 0;
  std::vector<struct MboAmplitude> mat(d * d);

  simpleTOpDenseMatrix(h, &sto, &mat[0]);
  for (int i = 0; i < d * d; ++i) {
    int row = i / d;
    int col = i % d;
    struct MboAmplitude expectedResult;
    if (row == col) {
      expectedResult.re = 1;
      expectedResult.im = 0;
    } else {
      expectedResult.re = 0;
      expectedResult.im = 0;
    }
    EXPECT_FLOAT_EQ(expectedResult.re, mat[i].re) << "(row, col == (" << row
                                                  << ", " << col << ")";
    EXPECT_FLOAT_EQ(expectedResult.im, mat[i].im) << "(row, col == (" << row
                                                  << ", " << col << ")";
  }

  destroySimpleTOp(&sto);
  mboProdSpaceDestroy(&h);
}

TEST(SimpleTOp, DenseMatrixIdentityTwoSpaces) {
  MboLocInd d = 2;
  MboProdSpace h = mboProdSpaceCreate(d);
  mboProdSpaceMul(h, &h);
  struct SimpleTOp sto;
  sto.numFactors = 0;
  sto.embeddings = 0;
  MboGlobInd dim = mboProdSpaceDim(h);
  std::vector<struct MboAmplitude> mat(dim * dim);

  simpleTOpDenseMatrix(h, &sto, &mat[0]);
  for (int i = 0; i < dim * dim; ++i) {
    int row = i / dim;
    int col = i % dim;
    struct MboAmplitude expectedResult;
    if (row == col) {
      expectedResult.re = 1;
      expectedResult.im = 0;
    } else {
      expectedResult.re = 0;
      expectedResult.im = 0;
    }
    EXPECT_FLOAT_EQ(expectedResult.re, mat[i].re) << "(row, col == (" << row
                                                  << ", " << col << ")";
    EXPECT_FLOAT_EQ(expectedResult.im, mat[i].im) << "(row, col == (" << row
                                                  << ", " << col << ")";
  }

  destroySimpleTOp(&sto);
  mboProdSpaceDestroy(&h);
}
