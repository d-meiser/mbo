#include <gtest/gtest.h>
#include <Tile.h>

TEST(Tile, SameTiles) {
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

TEST(Tile, DifferentTiles) {
  struct Tile t1 = {0, 1, 2, 3};
  struct Tile t2 = {1, 10, 0, 4};
  struct Tile intersection = tileIntersection(t1, t2);
  EXPECT_EQ(1, intersection.rmin);
  EXPECT_EQ(1, intersection.rmax);
  EXPECT_EQ(2, intersection.cmin);
  EXPECT_EQ(3, intersection.cmax);
}

TEST(Tile, IsEmptyFalse) {
  struct Tile t1 = {0, 2, 2, 3};
  struct Tile t2 = {1, 10, 0, 4};
  struct Tile intersection = tileIntersection(t1, t2);
  EXPECT_FALSE(tileIsEmpty(&intersection));
}

TEST(Tile, IsEmptyTrueR) {
  struct Tile t1 = {0, 1, 2, 3};
  struct Tile t2 = {1, 10, 0, 4};
  struct Tile intersection = tileIntersection(t1, t2);
  EXPECT_TRUE(tileIsEmpty(&intersection));
}

TEST(Tile, IsEmptyTrueC) {
  struct Tile t1 = {0, 2, 2, 3};
  struct Tile t2 = {1, 10, 0, 2};
  struct Tile intersection = tileIntersection(t1, t2);
  EXPECT_TRUE(tileIsEmpty(&intersection));
}

TEST(Tile, Divide) {
  struct Tile t1 = {1, 10, 0, 2};
  struct Tile t2 = {0, 2, 2, 3};
  struct Tile quotient = tileDivide(t1, t2);
  EXPECT_EQ(0, quotient.rmin);
  EXPECT_EQ(5, quotient.rmax);
  EXPECT_EQ(-2, quotient.cmin);
  EXPECT_EQ(0, quotient.cmax);
}

TEST(Tile, NumTilesContained) {
  struct Tile t1 = {1, 10, 0, 5};
  struct Tile t2 = {0, 2, 2, 3};
  MboGlobInd numTiles = numTilesContained(t1, t2);
  EXPECT_EQ(5, numTiles);
}

TEST(Tile, Advance) {
  struct Tile t = {0, 2, 2, 3};
  tileAdvance(3, &t);
  EXPECT_EQ(6, t.rmin);
  EXPECT_EQ(8, t.rmax);
  EXPECT_EQ(5, t.cmin);
  EXPECT_EQ(6, t.cmax);
}
