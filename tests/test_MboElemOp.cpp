#include <gtest/gtest.h>
#include <MboElemOp.h>
#include <MboAmplitude.h>

TEST(MboElemOp, Create) {
  MboElemOp a;
  mboElemOpCreate(&a);
  EXPECT_NE(a, (MboElemOp)0);
  mboElemOpDestroy(&a);
}

TEST(MboElemOp, AddTo) {
  MboElemOp a;
  mboElemOpCreate(&a);
  struct MboAmplitude alpha;
  alpha.re = 3.0;
  alpha.im = -5.0;
  mboElemOpAddTo(1, 2, &alpha, &a);
  mboElemOpDestroy(&a);
}

TEST(MboElemOp, Scale) {
  MboElemOp a;
  mboElemOpCreate(&a);
  struct MboAmplitude alpha, beta;
  alpha.re = 3.0;
  alpha.im = 1.0;
  mboElemOpAddTo(1, 2, &alpha, &a);
  beta.re = -1.0;
  beta.im = 8.0;
  mboElemOpScale(&beta, a);
  mboElemOpDestroy(&a);
}

TEST(MboElemOp, Plus) {
  MboElemOp a;
  mboElemOpCreate(&a);
  MboElemOp b = 0;
  mboElemOpCreate(&b);
  struct MboAmplitude alpha;
  alpha.re = 11;
  alpha.im = 22;
  mboElemOpAddTo(0, 1, &alpha, &a);
  mboElemOpPlus(a, &b);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &b);
  mboElemOpPlus(a, &b);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &a);
  mboElemOpAddTo(0, 1, &alpha, &b);
  mboElemOpPlus(a, &b);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &a);
  struct MboAmplitude beta;
  beta.re = -2.0;
  beta.im = 100;
  mboElemOpAddTo(5, 11, &beta, &b);
  mboElemOpPlus(a, &b);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);
}

TEST(MboElemOp, Mul) {
  struct MboAmplitude alpha;
  alpha.re = 3.0;
  alpha.im = 2.0;

  /* non-zero a * zero b */
  MboElemOp a;
  mboElemOpCreate(&a);
  MboElemOp b = 0;
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &a);
  mboElemOpMul(a, &b);
  EXPECT_EQ(mboElemOpNumEntries(a), 1);
  EXPECT_EQ(mboElemOpNumEntries(b), 0);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  /* zero a * non-zero b */
  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &b);
  mboElemOpMul(a, &b);
  EXPECT_EQ(mboElemOpNumEntries(a), 0);
  EXPECT_EQ(mboElemOpNumEntries(b), 0);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  /* non-zero a * non-zero b with zero product */
  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &a);
  mboElemOpAddTo(0, 1, &alpha, &b);
  mboElemOpMul(a, &b);
  EXPECT_EQ(mboElemOpNumEntries(a), 1);
  EXPECT_EQ(mboElemOpNumEntries(b), 0);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);

  /* non-zero a * non-zero b with non-zero product */
  mboElemOpCreate(&a);
  mboElemOpCreate(&b);
  mboElemOpAddTo(0, 1, &alpha, &a);
  mboElemOpAddTo(1, 0, &alpha, &b);
  mboElemOpMul(a, &b);
  EXPECT_EQ(mboElemOpNumEntries(a), 1);
  EXPECT_EQ(mboElemOpNumEntries(b), 1);
  mboElemOpDestroy(&a);
  mboElemOpDestroy(&b);
}

TEST(MboElemOp, Check) {
  MboElemOp o;

  mboElemOpCreate(&o);
  mboElemOpCheck(o);
  mboElemOpDestroy(&o);

  o = mboSigmaZ();
  mboElemOpCheck(o);
  mboElemOpDestroy(&o);
}

TEST(MboElemOp, SigmaPlus) {
  MboElemOp sp = mboSigmaPlus();
  EXPECT_EQ(mboElemOpNumEntries(sp), 1);
  mboElemOpDestroy(&sp);
}

TEST(MboElemOp, SigmaMinus) {
  MboElemOp sm = mboSigmaMinus();
  EXPECT_EQ(mboElemOpNumEntries(sm), 1);
  mboElemOpDestroy(&sm);
}

TEST(MboElemOp, SigmaZ) {
  MboElemOp sz = mboSigmaZ();
  EXPECT_EQ(mboElemOpNumEntries(sz), 2);
  mboElemOpDestroy(&sz);
}

TEST(MboElemOp, Eye) {
  MboElemOp e = mboEye(5);
  EXPECT_EQ(mboElemOpNumEntries(e), 5);
  mboElemOpDestroy(&e);
}

TEST(MboElemOp, NumOp) {
  MboElemOp num = mboNumOp(5);
  EXPECT_EQ(mboElemOpNumEntries(num), 4);
  mboElemOpDestroy(&num);
}

TEST(MboElemOp, AnnihilationOp) {
  MboElemOp a = mboAnnihilationOp(5);
  EXPECT_EQ(mboElemOpNumEntries(a), 4);
  mboElemOpDestroy(&a);
}

TEST(MboElemOp, CreationOp) {
  MboElemOp ad = mboCreationOp(5);
  EXPECT_EQ(mboElemOpNumEntries(ad), 4);
  mboElemOpDestroy(&ad);
}

