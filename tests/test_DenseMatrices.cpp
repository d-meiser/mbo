#include <gtest/gtest.h>
#include <MboTensorOp.h>
#include <MboNumOp.h>
#include <MboProdSpace.h>
#include <MboAmplitude.h>
#include <MboVec.h>

#include <iostream>

class TOpBuilder {
  public:
    virtual ~TOpBuilder() {}
    virtual MboNumOp build() const = 0;
    virtual TOpBuilder* copy() const = 0;
};

// Test Operators

class Null : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  Null* copy() const {
    return new Null(*this);
  }
};

class Identity : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpIdentity(h, &a);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  Identity* copy() const {
    return new Identity(*this);
  }
};

class SigmaPlus : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaPlus();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaPlus* copy() const {
    return new SigmaPlus(*this);
  }
};

class SigmaMinus : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaMinus();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaMinus* copy() const {
    return new SigmaMinus(*this);
  }
};

class SigmaZ : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaZ* copy() const {
    return new SigmaZ(*this);
  }
};

class AnnihilationOp : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(4);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboAnnihilationOp(4);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  AnnihilationOp* copy() const {
    return new AnnihilationOp(*this);
  }
};

class CreationOp : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(3);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboCreationOp(3);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  CreationOp* copy() const {
    return new CreationOp(*this);
  }
};

class NumOp : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboNumOp(2);
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  NumOp* copy() const {
    return new NumOp(*this);
  }
};

class SigmaPEnd : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaPlus();
    mboTensorOpAddTo(sp, 1, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaPEnd* copy() const {
    return new SigmaPEnd(*this);
  }
};

class SigmaZEnd : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 1, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaZEnd* copy() const {
    return new SigmaZEnd(*this);
  }
};

class SigmaZBegin : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 0, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
  }
  SigmaZBegin* copy() const {
    return new SigmaZBegin(*this);
  }
};

class SigmaZMiddle : public TOpBuilder {
 public:
  virtual MboNumOp build() const {
    MboProdSpace h = mboProdSpaceCreate(2);
    mboProdSpaceMul(h, &h);
    mboProdSpaceMul(h, &h);
    MboTensorOp a;
    mboTensorOpNull(h, &a);
    MboElemOp sp = mboSigmaZ();
    mboTensorOpAddTo(sp, 2, a);
    mboElemOpDestroy(&sp);
    mboProdSpaceDestroy(&h);
    MboNumOp ac;
    mboNumOpCompile(a, &ac);
    mboTensorOpDestroy(&a);
    return ac;
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
   MboNumOp build() const { return builder->build(); }

  private:
   BuilderWrapper();
   const TOpBuilder* builder;
};

class MboNumOpDenseMatrix
    : public ::testing::TestWithParam<BuilderWrapper> {};

// This function calculates the reference solution
static void computeMatrix(MboNumOp op, struct MboAmplitude* mat) {
  struct MboAmplitude one;
  one.re = 1.0;
  one.im = 0.0;
  struct MboAmplitude zero;
  zero.re = 0.0;
  zero.im = 0.0;
  MboGlobInd dim = mboProdSpaceDim(mboNumOpGetSpace(op));
  MboVec x;
  mboVecCreate(dim, &x);
  MboVec y;
  mboVecCreate(dim, &y);
  MboVec result;
  mboVecCreate(dim, &result);
  for (MboGlobInd i = 0; i < dim; ++i) {
    mboVecUnitVector(i, x);
    for (MboGlobInd j = 0; j < dim; ++j) {
      mboVecUnitVector(j, y);
      struct MboAmplitude *yarr;
      mboVecGetViewR(y, &yarr);
      struct MboAmplitude *resultarr;
      mboVecGetViewRW(result, &resultarr);
      mboNumOpMatVec(one, op, yarr, zero, resultarr);
      mboVecReleaseView(y, &yarr);
      mboVecReleaseView(result, &resultarr);
      mboVecDot(x, result, &mat[i * dim + j]);
    }
  }
  mboVecDestroy(&x);
  mboVecDestroy(&y);
  mboVecDestroy(&result);
}

static void printMatrix(struct MboAmplitude* mat, MboGlobInd dim,
                        const char* name) {
  std::cout << name << std::endl;
  std::cout << std::setprecision(1) << std::fixed;
  for (MboGlobInd i = 0; i < dim; ++i) {
    for (MboGlobInd j = 0; j < dim; ++j) {
      std::cout << mat[i * dim + j].re << " " << mat[i * dim + j].im
                << "  ";
    }
    std::cout << "\n";
  }
}

// The tests compare our results agains the reference solution
TEST_P(MboNumOpDenseMatrix, compAgainstNaive) {
  MboNumOp op = GetParam().build();
  MboGlobInd dim = mboProdSpaceDim(mboNumOpGetSpace(op));
  std::vector<struct MboAmplitude> mat(dim * dim);
  mboNumOpDenseMatrix(op, &mat[0]);
  std::vector<struct MboAmplitude> expected(dim * dim);
  computeMatrix(op, &expected[0]);
  for (MboGlobInd i = 0; i < dim; ++i) {
    for (MboGlobInd j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(expected[i * dim + j].re, mat[i * dim + j].re)
          << "Matrices differ in M(" << i << ", " << j << ").re";
      EXPECT_DOUBLE_EQ(expected[i * dim + j].im, mat[i * dim + j].im)
          << "Matrices differ in M(" << i << ", " << j << ").im";
    }
  }
  printMatrix(&expected[0], dim, "Expected matrix:");
  printMatrix(&mat[0], dim, "Actual matrix:");
  mboNumOpDestroy(&op);
}

INSTANTIATE_TEST_CASE_P(
    DenseMatrixTests, MboNumOpDenseMatrix,
    ::testing::Values(BuilderWrapper(Null()), BuilderWrapper(Identity()),
                      BuilderWrapper(SigmaPlus()), BuilderWrapper(SigmaMinus()),
                      BuilderWrapper(SigmaZ()),
                      BuilderWrapper(AnnihilationOp()),
                      BuilderWrapper(CreationOp()), BuilderWrapper(NumOp()),
                      BuilderWrapper(SigmaZEnd()), BuilderWrapper(SigmaPEnd()),
                      BuilderWrapper(SigmaZBegin()),
                      BuilderWrapper(SigmaZMiddle())));
