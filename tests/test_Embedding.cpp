#include <gtest/gtest.h>
#include <Embedding.h>

TEST(Embedding, Find) {
  int n = 5;
  struct Embedding *embeddings =
      (struct Embedding *)malloc(n * sizeof(*embeddings));

  embeddings[0].i = 3;
  embeddings[1].i = 4;
  embeddings[2].i = -1;
  embeddings[3].i = 3;
  embeddings[4].i = 0;

  EXPECT_EQ(findEmbedding(0, n, embeddings), 4);
  EXPECT_EQ(findEmbedding(3, n, embeddings), 0);
  EXPECT_EQ(findEmbedding(3, n - 1, embeddings + 1), 2);
  EXPECT_TRUE(findEmbedding(2, n, embeddings) < 0);
  EXPECT_EQ(findEmbedding(-1, n, embeddings), 2);

  free(embeddings);
}

TEST(Embedding, GatherIthEmbedding) {
  int i, n = 5;
  struct Embedding *embeddings =
      (struct Embedding *)malloc(n * sizeof(*embeddings));

  embeddings[0].i = 3;
  mboElemOpCreate(&embeddings[0].op);
  embeddings[1].i = 4;
  mboElemOpCreate(&embeddings[1].op);
  embeddings[2].i = 0;
  mboElemOpCreate(&embeddings[2].op);
  embeddings[3].i = 3;
  mboElemOpCreate(&embeddings[3].op);
  embeddings[4].i = 0;
  mboElemOpCreate(&embeddings[4].op);

  i = gatherIthEmbedding(0, &n, &embeddings);
  EXPECT_EQ(i, 2);
  EXPECT_EQ(n, 4);

  i = gatherIthEmbedding(4, &n, &embeddings);
  EXPECT_EQ(i, 1);
  EXPECT_EQ(n, 4);

  i = gatherIthEmbedding(3, &n, &embeddings);
  EXPECT_EQ(i, 0);
  EXPECT_EQ(n, 3);

  for (i = 0; i < n; ++i) {
    mboElemOpDestroy(&embeddings[i].op);
  }

  free(embeddings);
}

static int intCmp(const void *p1, const void *p2) {
  int *i1 = (int *)p1, *i2 = (int *)p2;
  return *i1 - *i2;
}

TEST(Embedding, GatherAllEmbeddings) {
  int n = 12, is[] = {2, 4, 3, 4, 4, 5, 6, 7, 9, 5, 6, 2}, isAfterGather[7];
  struct Embedding *embeddings =
      (struct Embedding *)malloc(n * sizeof(*embeddings));

  for (int i = 0; i < n; ++i) {
    embeddings[i].i = is[i];
    mboElemOpCreate(&embeddings[i].op);
  }

  gatherAllEmbeddings(&n, &embeddings);
  EXPECT_EQ(n, 7);
  for (int i = 0; i < n; ++i) {
    isAfterGather[i] = embeddings[i].i;
  }
  qsort(isAfterGather, n, sizeof(*isAfterGather), intCmp);
  EXPECT_EQ(isAfterGather[0], 2);
  EXPECT_EQ(isAfterGather[1], 3);
  EXPECT_EQ(isAfterGather[2], 4);
  EXPECT_EQ(isAfterGather[3], 5);
  EXPECT_EQ(isAfterGather[4], 6);
  EXPECT_EQ(isAfterGather[5], 7);
  EXPECT_EQ(isAfterGather[6], 9);

  for (int i = 0; i < n; ++i) {
    mboElemOpDestroy(&embeddings[i].op);
  }
  free(embeddings);
}

class ApplyEmbedding : public ::testing::Test {
  public:
  MboLocInd dims[4];
  MboGlobInd dim;
  struct MboAmplitude alpha, expectedResult, *x, *y;
  struct Embedding *embeddings;
  MboElemOp sz, sp;
  struct Tile tile;

  void SetUp() {
    dims[0] = 2;
    dims[1] = 3;
    dims[2] = 5;
    dims[3] = 2;
    x = (struct MboAmplitude *)malloc(computeBlockSize(4, dims) * sizeof(*x));
    y = (struct MboAmplitude *)malloc(computeBlockSize(4, dims) * sizeof(*y));
    sz = mboSigmaZ();
    sp = mboSigmaPlus();
    memset(x, 0, sizeof(*x) * computeBlockSize(4, dims));
    memset(y, 0, sizeof(*y) * computeBlockSize(4, dims));
    dim = computeBlockSize(4, dims);
  }

  void TearDown() {
    mboElemOpDestroy(&sz);
    mboElemOpDestroy(&sp);
    free(x);
    free(y);
  }
};

TEST_F(ApplyEmbedding, Identity) {
  tile.rmin = 0;
  tile.rmax = dim;
  tile.cmin = 0;
  tile.cmax = dim;
  applyEmbeddingsMask(0, 4, dims, dim, x[0], 0, 0, x, y, tile, &tile);
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    EXPECT_DOUBLE_EQ(y[i].re, 0);
    EXPECT_DOUBLE_EQ(y[i].im, 0);
  }
}

TEST_F(ApplyEmbedding, SzLast) {
  embeddings = (struct Embedding *)malloc(sizeof(*embeddings));
  embeddings[0].i = 3;
  embeddings[0].op = sz;
  alpha.re = 2.0;
  alpha.im = 3.0;
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    x[i].re = 1.0;
    x[i].im = 0.0;
  }
  tile.rmin = 0;
  tile.rmax = dim;
  tile.cmin = 0;
  tile.cmax = dim;
  applyEmbeddingsMask(0, 4, dims, dim, alpha, 1, embeddings, x, y, tile, &tile);
  for (int i = 0; i < dim; ++i) {
    if (i & 1l) {
      expectedResult.re = alpha.re;
      expectedResult.im = alpha.im;
    } else {
      expectedResult.re = -alpha.re;
      expectedResult.im = -alpha.im;
    }
    EXPECT_DOUBLE_EQ(expectedResult.re, y[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(expectedResult.im, y[i].im) << " i == " << i;
  }
  free(embeddings);
}

TEST_F(ApplyEmbedding, SzInterior) {
  embeddings = (struct Embedding *)malloc(sizeof(*embeddings));
  embeddings[0].i = 1;
  embeddings[0].op = sz;
  alpha.re = 2.0;
  alpha.im = 3.0;
  for (int i = 0; i < dim; ++i) {
    x[i].re = 1.0;
    x[i].im = 0.0;
  }
  for (int i = 0; i < dim; ++i) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }
  tile.rmin = 0;
  tile.rmax = dim;
  tile.cmin = 0;
  tile.cmax = dim;
  applyEmbeddingsMask(0, 4, dims, dim, alpha, 1, embeddings, x, y, tile, &tile);
  for (int i = 0; i < dim; ++i) {
    switch ((i / computeBlockSize(2, dims + 2)) % dims[1]) {
      case 0:
        expectedResult.re = -alpha.re;
        expectedResult.im = -alpha.im;
        break;
      case 1:
        expectedResult.re = alpha.re;
        expectedResult.im = alpha.im;
        break;
      default:
        expectedResult.re = 0;
        expectedResult.im = 0;
    }
    EXPECT_DOUBLE_EQ(y[i].re, expectedResult.re);
    EXPECT_DOUBLE_EQ(y[i].im, expectedResult.im);
  }
  free(embeddings);
}

TEST_F(ApplyEmbedding, TwoFactors) {
  embeddings = (struct Embedding *)malloc(2 * sizeof(*embeddings));
  embeddings[0].i = 0;
  embeddings[0].op = sz;
  embeddings[1].i = 2;
  embeddings[1].op = sp;
  alpha.re = 2.0;
  alpha.im = 3.0;
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    x[i].re = 1.0;
    x[i].im = 0.0;
  }
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }
  applyEmbeddingsMask(0, 4, dims, computeBlockSize(4, dims), alpha, 2, embeddings,
                  x, y, tile, &tile);
  /* Make sure that dim of first space is 2.  Otherwise the expected
   * results are not calculated correctly. */
  ASSERT_EQ(dims[0], 2);
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    if ((i / computeBlockSize(1, dims + 3)) % dims[2] == 1) {
      if ((i / computeBlockSize(3, dims + 1)) % dims[0]) {
        expectedResult.re = alpha.re;
        expectedResult.im = alpha.im;
      } else {
        expectedResult.re = -alpha.re;
        expectedResult.im = -alpha.im;
        break;
      }
    } else {
      expectedResult.re = 0;
      expectedResult.im = 0;
    }
    EXPECT_DOUBLE_EQ(y[i].re, expectedResult.re);
    EXPECT_DOUBLE_EQ(y[i].im, expectedResult.im);
  }
  free(embeddings);
}

TEST(Embedding, SortEmbeddings) {
  int i;
  int numEmbeddings;
  struct Embedding *embeddings;

  numEmbeddings = 0;
  embeddings = 0;
  sortEmbeddings(numEmbeddings, embeddings);
  EXPECT_EQ(numEmbeddings, 0);

  numEmbeddings = 1;
  embeddings = (struct Embedding *)malloc(numEmbeddings * sizeof(*embeddings));
  embeddings[0].i = 3;
  for (i = 0; i < numEmbeddings; ++i) {
    mboElemOpCreate(&embeddings[i].op);
  }
  sortEmbeddings(numEmbeddings, embeddings);
  EXPECT_EQ(numEmbeddings, 1);
  EXPECT_EQ(embeddings[0].i, 3);
  for (i = 0; i < numEmbeddings; ++i) {
    destroyEmbedding(&embeddings[i]);
  }
  free(embeddings);

  /* With duplicates */
  numEmbeddings = 6;
  embeddings = (struct Embedding *)malloc(numEmbeddings * sizeof(*embeddings));
  embeddings[0].i = 3;
  embeddings[1].i = 4;
  embeddings[2].i = 2;
  embeddings[3].i = 4;
  embeddings[4].i = 4;
  embeddings[5].i = 0;
  for (i = 0; i < numEmbeddings; ++i) {
    mboElemOpCreate(&embeddings[i].op);
  }
  sortEmbeddings(numEmbeddings, embeddings);
  for (i = 1; i < numEmbeddings; ++i) {
    EXPECT_TRUE(embeddings[i].i >= embeddings[i - 1].i);
  }
  for (i = 0; i < numEmbeddings; ++i) {
    destroyEmbedding(&embeddings[i]);
  }
  free(embeddings);

  /* Without duplicates we end up with an ordered array */
  numEmbeddings = 6;
  embeddings = (struct Embedding *)malloc(numEmbeddings * sizeof(*embeddings));
  embeddings[0].i = 3;
  embeddings[1].i = 4;
  embeddings[2].i = 2;
  embeddings[3].i = 12;
  embeddings[4].i = 15;
  embeddings[5].i = 7;
  for (i = 0; i < numEmbeddings; ++i) {
    mboElemOpCreate(&embeddings[i].op);
  }
  sortEmbeddings(numEmbeddings, embeddings);
  for (i = 1; i < numEmbeddings; ++i) {
    EXPECT_TRUE(embeddings[i].i > embeddings[i - 1].i);
  }
  for (i = 0; i < numEmbeddings; ++i) {
    destroyEmbedding(&embeddings[i]);
  }
  free(embeddings);
}

TEST(Embedding, ApplyEmbeddingsRowRange) {
  MboLocInd dims[] = {2, 3};
  struct Tile tile;
  
  struct MboAmplitude tmp;
  tmp.re = 1;
  tmp.im = 0;
  std::vector<struct MboAmplitude> x(computeBlockSize(2, dims));
  std::fill(x.begin(), x.end(), tmp);
  std::vector<struct MboAmplitude> y(computeBlockSize(2, dims));
  tmp.re = 0;
  tmp.im = 0;
  std::fill(y.begin(), y.end(), tmp);

  MboElemOp sz = mboSigmaZ();
  std::vector<struct Embedding> embeddings(1);
  embeddings[0].i = 0;
  embeddings[0].op = sz;
  struct MboAmplitude alpha = {2.0, 3.0};
  tile.rmin = 0;
  tile.rmax = 3;
  tile.cmin = 0;
  tile.cmax = computeBlockSize(2, dims);
  applyEmbeddingsMask(0, 2, dims, computeBlockSize(2, dims), alpha,
      1, &embeddings[0], &x[0], &y[0], tile, &tile);
  EXPECT_DOUBLE_EQ(-alpha.re, y[0].re) << " i == " << 0;
  EXPECT_DOUBLE_EQ(-alpha.im, y[0].im) << " i == " << 0;
  EXPECT_DOUBLE_EQ(-alpha.re, y[1].re) << " i == " << 1;
  EXPECT_DOUBLE_EQ(-alpha.im, y[1].im) << " i == " << 1;
  EXPECT_DOUBLE_EQ(-alpha.re, y[2].re) << " i == " << 2;
  EXPECT_DOUBLE_EQ(-alpha.im, y[2].im) << " i == " << 2;
  for (MboLocInd i = 3; i < computeBlockSize(2, dims); ++i) {
    EXPECT_DOUBLE_EQ(0, y[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(0, y[i].im) << " i == " << i;
  }
  mboElemOpDestroy(&sz);
}

TEST(Embedding, ApplyEmbeddingsRmaxOutOfRange) {
  MboLocInd dims[] = {2, 3};
  struct Tile tile, fullTile;

  struct MboAmplitude tmp;
  tmp.re = 1;
  tmp.im = 0;
  std::vector<struct MboAmplitude> x(computeBlockSize(2, dims));
  std::fill(x.begin(), x.end(), tmp);
  std::vector<struct MboAmplitude> y(computeBlockSize(2, dims));
  tmp.re = 0;
  tmp.im = 0;
  std::fill(y.begin(), y.end(), tmp);

  MboElemOp sz = mboSigmaZ();
  std::vector<struct Embedding> embeddings(1);
  embeddings[0].i = 0;
  embeddings[0].op = sz;
  struct MboAmplitude alpha = {2.0, 3.0};
  tile.rmin = 4;
  tile.rmax = 7;
  tile.cmin = 0;
  tile.cmax = computeBlockSize(2, dims);
  fullTile.rmin = 0;
  fullTile.rmax = computeBlockSize(2, dims);
  fullTile.cmin = 0;
  fullTile.cmax = computeBlockSize(2, dims);
  applyEmbeddingsMask(0, 2, dims, computeBlockSize(2, dims), alpha,
      1, &embeddings[0], &x[0], &y[0] + 4, fullTile, &tile);
  EXPECT_DOUBLE_EQ(alpha.re, y[4].re) << " i == " << 4;
  EXPECT_DOUBLE_EQ(alpha.im, y[4].im) << " i == " << 4;
  EXPECT_DOUBLE_EQ(alpha.re, y[5].re) << " i == " << 5;
  EXPECT_DOUBLE_EQ(alpha.im, y[5].im) << " i == " << 5;
  for (MboLocInd i = 0; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(0, y[i].re) << " i == " << i;
    EXPECT_DOUBLE_EQ(0, y[i].im) << " i == " << i;
  }
  mboElemOpDestroy(&sz);
}

struct ApplyLeafMask : public ::testing::Test {
  struct Tile tile, mask;
  struct MboAmplitude alpha;
  std::vector<struct MboAmplitude> x;
  std::vector<struct MboAmplitude> y;

  void SetUp() {
    alpha.re = 1.0;
    alpha.im = 0.0;
    tile.rmin = 0;
    tile.cmin = 0;
    tile.rmax = 0;
    tile.cmax = 0;
    mask.rmin = 0;
    mask.cmin = 0;
    mask.rmax = 0;
    mask.cmax = 0;
  }
};

TEST_F(ApplyLeafMask, Empty) {
  y.resize(1);
  y[0].re = 0;
  y[0].im = 0;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(0, y[0].re);
}

TEST_F(ApplyLeafMask, TileEmpty) {
  y.resize(1);
  y[0].re = 0;
  y[0].im = 0;
  mask.rmax = 1;
  mask.cmax = 1;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(0, y[0].re);
}

TEST_F(ApplyLeafMask, MaskEmpty) {
  y.resize(1);
  y[0].re = 0;
  y[0].im = 0;
  tile.rmax = 1;
  tile.cmax = 1;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(0, y[0].re);
}

TEST_F(ApplyLeafMask, MaskAndTileEqual) {
  y.resize(1);
  y[0].re = 12.0;
  y[0].im = 20.0;
  x = y;
  tile.rmax = 1;
  tile.cmax = 1;
  mask = tile;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(24.0, y[0].re);
  EXPECT_DOUBLE_EQ(40.0, y[0].im);
}

TEST_F(ApplyLeafMask, TileFullyInsideMask) {
  y.resize(4);
  struct MboAmplitude yIn;
  yIn.re = 12.0;
  yIn.im = 20.0;
  std::fill(y.begin(), y.end(), yIn);
  struct MboAmplitude xIn;
  xIn.re = 2.0;
  xIn.im = 5.0;
  x.resize(3);
  std::fill(x.begin(), x.end(), xIn);
  tile.rmin = 1;
  tile.cmin = 1;
  tile.rmax = 2;
  tile.cmax = 2;
  mask.rmax = 4;
  mask.cmax = 3;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(yIn.re, y[0].re);
  EXPECT_DOUBLE_EQ(yIn.im, y[0].im);
  EXPECT_DOUBLE_EQ(yIn.re + xIn.re, y[1].re);
  EXPECT_DOUBLE_EQ(yIn.im + xIn.im, y[1].im);
  EXPECT_DOUBLE_EQ(yIn.re, y[2].re);
  EXPECT_DOUBLE_EQ(yIn.im, y[2].im);
  EXPECT_DOUBLE_EQ(yIn.re, y[3].re);
  EXPECT_DOUBLE_EQ(yIn.im, y[3].im);
}

TEST_F(ApplyLeafMask, MaskFullyInsideTile) {
  y.resize(1);
  struct MboAmplitude yIn;
  yIn.re = 12.0;
  yIn.im = 20.0;
  std::fill(y.begin(), y.end(), yIn);
  struct MboAmplitude xIn;
  xIn.re = 2.0;
  xIn.im = 5.0;
  x.resize(3);
  std::fill(x.begin(), x.end(), xIn);
  tile.rmax = 4;
  tile.cmax = 5;
  mask.rmin = 1;
  mask.cmin = 0;
  mask.rmax = 2;
  mask.cmax = 3;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  EXPECT_DOUBLE_EQ(yIn.re + xIn.re, y[0].re);
  EXPECT_DOUBLE_EQ(yIn.im + xIn.im, y[0].im);
}

TEST_F(ApplyLeafMask, MaskBeyondRowRange) {
  y.resize(4);
  struct MboAmplitude yIn;
  yIn.re = 12.0;
  yIn.im = 20.0;
  std::fill(y.begin(), y.end(), yIn);
  struct MboAmplitude xIn;
  xIn.re = 2.0;
  xIn.im = 5.0;
  x.resize(3);
  std::fill(x.begin(), x.end(), xIn);
  tile.rmax = 4;
  tile.cmax = 5;
  mask.rmin = 1;
  mask.cmin = 0;
  mask.rmax = 5;
  mask.cmax = 3;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  for (int i = 0; i < 2; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re + xIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im + xIn.im, y[i].im) << "i == " << i;
  }
  for (int i = 2; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im, y[i].im) << "i == " << i;
  }
}

TEST_F(ApplyLeafMask, MaskBeyondRowAndColRange) {
  y.resize(3);
  struct MboAmplitude yIn;
  yIn.re = 12.0;
  yIn.im = 20.0;
  std::fill(y.begin(), y.end(), yIn);
  struct MboAmplitude xIn;
  xIn.re = 2.0;
  xIn.im = 5.0;
  x.resize(6);
  std::fill(x.begin(), x.end(), xIn);
  tile.rmin = 3;
  tile.cmin = 3;
  tile.rmax = 6;
  tile.cmax = 6;
  mask.rmin = 4;
  mask.cmin = 0;
  mask.rmax = 7;
  mask.cmax = 6;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  for (int i = 0; i < 2; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re + xIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im + xIn.im, y[i].im) << "i == " << i;
  }
  for (int i = 2; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im, y[i].im) << "i == " << i;
  }
}

TEST_F(ApplyLeafMask, MaskBeyondRowAndColRangeTranspose) {
  y.resize(6);
  struct MboAmplitude yIn;
  yIn.re = 12.0;
  yIn.im = 20.0;
  std::fill(y.begin(), y.end(), yIn);
  struct MboAmplitude xIn;
  xIn.re = 2.0;
  xIn.im = 5.0;
  x.resize(3);
  std::fill(x.begin(), x.end(), xIn);
  tile.rmin = 3;
  tile.cmin = 3;
  tile.rmax = 6;
  tile.cmax = 6;
  mask.rmin = 0;
  mask.cmin = 4;
  mask.rmax = 6;
  mask.cmax = 7;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  for (int i = 0; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im, y[i].im) << "i == " << i;
  }
  for (int i = 4; i < 6; ++i) {
    EXPECT_DOUBLE_EQ(yIn.re + xIn.re, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(yIn.im + xIn.im, y[i].im) << "i == " << i;
  }
}

TEST_F(ApplyLeafMask, MaskRowBand) {
  y.resize(3);
  struct MboAmplitude yIn;
  yIn.re = 0.0;
  yIn.im = 0.0;
  std::fill(y.begin(), y.end(), yIn);
  x.resize(6);
  for (int i = 0; i < 6; ++i) {
    x[i].re = i;
    x[i].im = 0;
  }
  tile.rmin = 0;
  tile.cmin = 3;
  tile.rmax = 3;
  tile.cmax = 6;
  mask.rmin = 0;
  mask.cmin = 0;
  mask.rmax = 3;
  mask.cmax = 6;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(3.0 + i, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(0.0, y[i].im) << "i == " << i;
  }
}

TEST_F(ApplyLeafMask, MaskRowBand2) {
  y.resize(2);
  struct MboAmplitude yIn;
  yIn.re = 0.0;
  yIn.im = 0.0;
  std::fill(y.begin(), y.end(), yIn);
  x.resize(6);
  for (int i = 0; i < 6; ++i) {
    x[i].re = i;
    x[i].im = 0;
  }
  tile.rmin = 0;
  tile.cmin = 3;
  tile.rmax = 3;
  tile.cmax = 6;
  mask.rmin = 1;
  mask.cmin = 0;
  mask.rmax = 3;
  mask.cmax = 6;
  applyLeafMask(alpha, &x[0], &y[0], &tile, &mask);
  for (int i = 0; i < 2; ++i) {
    EXPECT_DOUBLE_EQ(4.0 + i, y[i].re) << "i == " << i;
    EXPECT_DOUBLE_EQ(0.0, y[i].im) << "i == " << i;
  }
}
