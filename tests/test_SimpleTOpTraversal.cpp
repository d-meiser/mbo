#include <gtest/gtest.h>
#include <SimpleTOpTraversal.h>

TEST(Tile, SameTiles)
{
	struct Tile t;
	t.rmin = 0;
	t.rmax = 1;
	t.cmin = 2;
	t.cmax = 3;
	struct Tile intersection = tileIntersection(t, t);
	EXPECT_EQ(t.rmin, intersection.rmin);
	EXPECT_EQ(t.rmax, intersection.rmax);
	EXPECT_EQ(t.cmin, intersection.cmin);
	EXPECT_EQ(t.cmax, intersection.cmax);
}

TEST(Tile, DifferentTiles)
{
	struct Tile t1 = {0, 1, 2, 3};
	struct Tile t2 = {1, 10, 0, 4};
        struct Tile intersection = tileIntersection(t1, t2);
        EXPECT_EQ(1, intersection.rmin);
	EXPECT_EQ(1, intersection.rmax);
	EXPECT_EQ(2, intersection.cmin);
	EXPECT_EQ(3, intersection.cmax);
}

TEST(Tile, IsEmptyFalse)
{
	struct Tile t1 = {0, 2, 2, 3};
	struct Tile t2 = {1, 10, 0, 4};
        struct Tile intersection = tileIntersection(t1, t2);
        EXPECT_FALSE(tileIsEmpty(&intersection));
}

TEST(Tile, IsEmptyTrueR)
{
	struct Tile t1 = {0, 1, 2, 3};
	struct Tile t2 = {1, 10, 0, 4};
        struct Tile intersection = tileIntersection(t1, t2);
        EXPECT_TRUE(tileIsEmpty(&intersection));
}

TEST(Tile, IsEmptyTrueC)
{
	struct Tile t1 = {0, 2, 2, 3};
	struct Tile t2 = {1, 10, 0, 2};
        struct Tile intersection = tileIntersection(t1, t2);
        EXPECT_TRUE(tileIsEmpty(&intersection));
}
