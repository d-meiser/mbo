#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboProdSpace.h>
#include <MboAmplitude.h>
#include <MboVec.h>

class TOpBuilder {
  public:
    virtual ~TOpBuilder() {}
    virtual MboTensorOp build() const = 0;
    virtual TOpBuilder* copy() const = 0;
};

class NullBuilder : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    mboProdSpaceDestroy(&h);
    return a;
  }
  NullBuilder* copy() const {
    return new NullBuilder(*this);
  }
};

class IdentityBuilder : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpIdentity(h, &a);
    mboProdSpaceDestroy(&h);
    return a;
  }
  IdentityBuilder* copy() const {
    return new IdentityBuilder(*this);
  }
};

class BuilderWrapper {
  public:
   BuilderWrapper(const TOpBuilder& b) : builder(b.copy()) {}
   BuilderWrapper(const BuilderWrapper& other)
       : builder(other.builder->copy()) {}
   BuilderWrapper& operator=(const BuilderWrapper& other) {
     if (this != &other) {
       delete this->builder;
       this->builder = other.builder->copy();
     }
     return *this;
   }
   ~BuilderWrapper() { delete builder; }
   MboTensorOp build() const { return builder->build(); }

  private:
   BuilderWrapper();
   const TOpBuilder* builder;
};

class MboTensorOpDenseMatrix
    : public ::testing::TestWithParam<BuilderWrapper> {};

static void computeMatrix(MboTensorOp op, struct MboAmplitude* mat) {
  struct MboAmplitude one;
  one.re = 1.0;
  one.im = 0.0;
  struct MboAmplitude zero;
  zero.re = 0.0;
  zero.im = 0.0;
  MboGlobInd dim = mboProdSpaceDim(mboTensorOpGetSpace(op));
  MboVec x;
  mboVecCreate(dim, &x);
  MboVec y;
  mboVecCreate(dim, &y);
  for (MboGlobInd i = 0; i < dim; ++i) {
    mboVecUnitVector(i, x);
    for (MboGlobInd j = 0; j < dim; ++j) {
      mboVecUnitVector(j, y);
      mboTensorOpMatVec(&one, op, y, &zero, y);
      mboVecDot(x, y, &mat[i * dim + j]);
    }
  }
  mboVecDestroy(&x);
  mboVecDestroy(&y);
}

TEST_P(MboTensorOpDenseMatrix, compAgainstNaive) {
  MboTensorOp op = GetParam().build();
  MboGlobInd dim = mboProdSpaceDim(mboTensorOpGetSpace(op));
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboTensorOpDenseMatrix(op, &mat[0]);
  std::vector<struct MboAmplitude> expected(dim * dim);
  computeMatrix(op, &expected[0]);
  for (MboGlobInd i = 0; i < dim; ++i) {
    for (MboGlobInd j = 0; j < dim; ++j) {
      EXPECT_FLOAT_EQ(mat[i * dim + j].re, expected[i * dim + j].re)
          << "Matrices differ in M(" << i << ", " << j << ").re";
      EXPECT_FLOAT_EQ(mat[i * dim + j].im, expected[i * dim + j].im)
          << "Matrices differ in M(" << i << ", " << j << ").im";
    }
  }
  mboTensorOpDestroy(&op);
}

INSTANTIATE_TEST_CASE_P(DenseMatrixTests, MboTensorOpDenseMatrix,
                        ::testing::Values(BuilderWrapper(NullBuilder()),
                                          BuilderWrapper(IdentityBuilder())));
