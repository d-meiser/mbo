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

// Test Operators

class Null : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    mboProdSpaceDestroy(&h);
    return a;
  }
  Null* copy() const {
    return new Null(*this);
  }
};

class Identity : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpIdentity(h, &a);
    mboProdSpaceDestroy(&h);
    return a;
  }
  Identity* copy() const {
    return new Identity(*this);
  }
};

class SigmaPlus : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaPlus();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaPlus* copy() const {
    return new SigmaPlus(*this);
  }
};

class SigmaMinus : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaMinus();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaMinus* copy() const {
    return new SigmaMinus(*this);
  }
};

class SigmaZ : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaZ* copy() const {
    return new SigmaZ(*this);
  }
};

class AnnihilationOp : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(4);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboAnnihilationOp(4);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  AnnihilationOp* copy() const {
    return new AnnihilationOp(*this);
  }
};

class CreationOp : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(3);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboCreationOp(3);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  CreationOp* copy() const {
    return new CreationOp(*this);
  }
};

class NumOp : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboNumOp(2);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  NumOp* copy() const {
    return new NumOp(*this);
  }
};

class SigmaPEnd : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaPlus();
    mboTensorOpAddTo(sp, 1, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaPEnd* copy() const {
    return new SigmaPEnd(*this);
  }
};

class SigmaZEnd : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 1, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaZEnd* copy() const {
    return new SigmaZEnd(*this);
  }
};

class SigmaZBegin : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaZBegin* copy() const {
    return new SigmaZBegin(*this);
  }
};

class SigmaZMiddle : public TOpBuilder {
 public:
  virtual MboTensorOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 2, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    return a;
  }
  SigmaZMiddle* copy() const {
    return new SigmaZMiddle(*this);
  }
};

// Some boiler plate for feeding the test operators into the test cases.

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

// This function calculates the reference solution
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

// The tests compare our results agains the reference solution
TEST_P(MboTensorOpDenseMatrix, compAgainstNaive) {
  MboTensorOp op = GetParam().build();
  MboGlobInd dim = mboProdSpaceDim(mboTensorOpGetSpace(op));
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboTensorOpDenseMatrix(op, &mat[0]);
  std::vector<struct MboAmplitude> expected(dim * dim);
  computeMatrix(op, &expected[0]);
  for (MboGlobInd i = 0; i < dim; ++i) {
    for (MboGlobInd j = 0; j < dim; ++j) {
      EXPECT_FLOAT_EQ(expected[i * dim + j].re, mat[i * dim + j].re)
          << "Matrices differ in M(" << i << ", " << j << ").re";
      EXPECT_FLOAT_EQ(expected[i * dim + j].im, mat[i * dim + j].im)
          << "Matrices differ in M(" << i << ", " << j << ").im";
    }
  }
  mboTensorOpDestroy(&op);
}

INSTANTIATE_TEST_CASE_P(DenseMatrixTests, MboTensorOpDenseMatrix,
                        ::testing::Values(BuilderWrapper(Null()),
                                          BuilderWrapper(Identity()),
                                          BuilderWrapper(SigmaPlus()),
                                          BuilderWrapper(SigmaMinus()),
                                          BuilderWrapper(SigmaZ()),
                                          BuilderWrapper(AnnihilationOp()),
                                          BuilderWrapper(CreationOp()),
                                          BuilderWrapper(NumOp()),
                                          BuilderWrapper(SigmaZEnd()),
                                          BuilderWrapper(SigmaPEnd()),
					  BuilderWrapper(SigmaZBegin()),
					  BuilderWrapper(SigmaZMiddle())));
