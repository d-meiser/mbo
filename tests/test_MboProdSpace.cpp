#include <gtest/gtest.h>
#include <algorithm>
#include <MboProdSpace.h>

TEST(MboProdSpace, InternalTests) { EXPECT_EQ(0, mboProdSpaceTest()); }

TEST(MboProdSpace, mboProdSpaceCreate) {
  MboProdSpace sp;
  sp = mboProdSpaceCreate(2);
  EXPECT_NE(sp, (MboProdSpace)0);
  mboProdSpaceDestroy(&sp);

  sp = mboProdSpaceCreate(0);
  EXPECT_NE(sp, (MboProdSpace)0);
  mboProdSpaceDestroy(&sp);
}

TEST(MboProdSpace, mboProdSpaceDestroy) {
  MboProdSpace sp = mboProdSpaceCreate(1);
  mboProdSpaceDestroy(&sp);
  EXPECT_EQ(sp, (MboProdSpace)0);
}

class MboProdSpaceMulTest
    : public ::testing::TestWithParam<std::vector<MboLocInd> > {};

static MboGlobInd dim(const std::vector<MboLocInd>& dims) {
  MboGlobInd result = 1;
  for (size_t i = 0; i < dims.size(); ++i) {
    result *= dims[i] > 0 ? dims[i] : 1;
  }
  return result;
}

bool nonzero(MboLocInd i) { return i != 0; }

MboProdSpace buildProdSpace(const std::vector<MboLocInd>& dims) {
  MboProdSpace sp = mboProdSpaceCreate(0);
  for (size_t i = 0; i < dims.size(); ++i) {
    MboProdSpace tmp = mboProdSpaceCreate(dims[i]);
    mboProdSpaceMul(tmp, &sp);
    mboProdSpaceDestroy(&tmp);
  }
  return sp;
}

TEST_P(MboProdSpaceMulTest, multiplication) {
  std::vector<MboLocInd> dims = GetParam();
  MboProdSpace sp = buildProdSpace(dims);
  EXPECT_EQ(mboProdSpaceDim(sp), dim(dims));
  mboProdSpaceDestroy(&sp);
}

std::vector<MboLocInd> buildRandomDims(int n) {
  std::vector<MboLocInd> dims;
  for (int i = 0; i < n; ++i) {
    dims.push_back(rand() % 10);
  }
  return dims;
}
std::vector<std::vector<MboLocInd> > testDims() {
  std::vector<std::vector<MboLocInd> > cases;
  cases.push_back(std::vector<MboLocInd>());
  cases.push_back(std::vector<MboLocInd>(1));
  cases.push_back(std::vector<MboLocInd>(3));
  cases.push_back(std::vector<MboLocInd>(1, 3));
  cases.push_back(std::vector<MboLocInd>(2, 3));
  cases.push_back(buildRandomDims(5));
  return cases;
}
INSTANTIATE_TEST_CASE_P(MulTests, MboProdSpaceMulTest,
                        ::testing::ValuesIn(testDims()));

std::vector<MboLocInd> manyDims(20, 2);
INSTANTIATE_TEST_CASE_P(BuildSpace, MboProdSpaceMulTest,
                        ::testing::Values(manyDims));

std::vector<MboLocInd> largeDims(2, 0x7FFFFFFF);
INSTANTIATE_TEST_CASE_P(LargeDims, MboProdSpaceMulTest,
                        ::testing::Values(largeDims));

TEST(MboProdSpace, MultiplyWithSelf)
{
  MboGlobInd d = 2;
  MboProdSpace h = mboProdSpaceCreate((MboLocInd)d);
  EXPECT_EQ(mboProdSpaceDim(h), d);
  for (int i = 0; i < 3; ++i) {
    d *= d;
    mboProdSpaceMul(h, &h);
    EXPECT_EQ(mboProdSpaceDim(h), d);
  }
  mboProdSpaceDestroy(&h);
}

TEST(MboProdSpace, GetDims)
{
  std::vector<MboLocInd> dims;
  dims.push_back(1);
  dims.push_back(5);
  dims.push_back(0);
  dims.push_back(3);
  dims.push_back(7);
  MboProdSpace h = buildProdSpace(dims);
  std::vector<MboLocInd> hDims(mboProdSpaceSize(h));
  EXPECT_EQ(4, hDims.size());
  mboProdSpaceGetDims(h, hDims.size(), hDims.empty() ? 0 : &hDims[0]);
  EXPECT_EQ(7, hDims[0]);
  EXPECT_EQ(3, hDims[1]);
  EXPECT_EQ(5, hDims[2]);
  EXPECT_EQ(1, hDims[3]);
  mboProdSpaceDestroy(&h);
}

TEST(MboProdSpace, CheckRand)
{
  MboProdSpace h = buildProdSpace(buildRandomDims(3));
  EXPECT_EQ(mboProdSpaceCheck(h), 0);
  mboProdSpaceDestroy(&h);
}

TEST(MboProdSpace, CheckEmpty)
{
  MboProdSpace h = buildProdSpace(std::vector<MboLocInd>());
  EXPECT_EQ(mboProdSpaceCheck(h), 0);
  mboProdSpaceDestroy(&h);
}

TEST(MboProdSpace, CheckFail)
{
  std::vector<MboLocInd> dims = buildRandomDims(3);
  dims[1] = -3;
  MboProdSpace h = buildProdSpace(dims);
  EXPECT_NE(mboProdSpaceCheck(h), 0);
  mboProdSpaceDestroy(&h);
}
