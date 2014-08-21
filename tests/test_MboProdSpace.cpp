#include <gtest/gtest.h>
#include <MboProdSpace.h>

TEST(MboProdSpace, InternalTests)
{
	EXPECT_EQ(0, mboProdSpaceTest());
}

TEST(MboProdSpace, mboCreateProdSpace)
{
	MboProdSpace sp;
	sp = mboProdSpaceCreate(2);
	EXPECT_NE(sp, (MboProdSpace)0);
	mboProdSpaceDestroy(&sp);

	sp = mboProdSpaceCreate(0);
	EXPECT_NE(sp, (MboProdSpace)0);
	mboProdSpaceDestroy(&sp);
}

TEST(MboProdSpace, mboCreateProdDestroy)
{
	MboProdSpace sp = mboProdSpaceCreate(1);
	mboProdSpaceDestroy(&sp);
	EXPECT_EQ(sp, (MboProdSpace)0);
}