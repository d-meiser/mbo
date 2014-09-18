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

TEST(EMBEDDING, GatherIthEmbedding) {
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

TEST(Embedding, ApplyEmbeddings) {
  MboLocInd dims[] = {2, 3, 5, 2};
  MboGlobInd blockSize;
  struct MboAmplitude alpha, expectedResult, *x, *y;
  struct Embedding *embeddings;
  MboElemOp sz, sp;

  x = (struct MboAmplitude *)malloc(computeBlockSize(4, dims) * sizeof(*x));
  y = (struct MboAmplitude *)malloc(computeBlockSize(4, dims) * sizeof(*y));
  sz = mboSigmaZ();
  sp = mboSigmaPlus();

  memset(x, 0, sizeof(*x) * computeBlockSize(4, dims));
  memset(y, 0, sizeof(*y) * computeBlockSize(4, dims));

  blockSize = computeBlockSize(3, dims + 1);
  applyEmbeddings(0, 4, dims, blockSize, x[0], 0, 0, x, y, 0,
                  computeBlockSize(4, dims));
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    EXPECT_FLOAT_EQ(y[i].re, 0);
    EXPECT_FLOAT_EQ(y[i].im, 0);
  }

  embeddings = (struct Embedding *)malloc(sizeof(*embeddings));
  embeddings[0].i = 3;
  embeddings[0].op = sz;
  alpha.re = 2.0;
  alpha.im = 3.0;
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    x[i].re = 1.0;
    x[i].im = 0.0;
  }
  applyEmbeddings(0, 4, dims, computeBlockSize(4, dims), alpha, 1, embeddings,
                  x, y, 0, computeBlockSize(4, dims));
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
    if (i & 1l) {
      expectedResult.re = alpha.re;
      expectedResult.im = alpha.im;
    } else {
      expectedResult.re = -alpha.re;
      expectedResult.im = -alpha.im;
    }
    EXPECT_FLOAT_EQ(y[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(y[i].im, expectedResult.im);
  }
  free(embeddings);

  embeddings = (struct Embedding *)malloc(sizeof(*embeddings));
  embeddings[0].i = 1;
  embeddings[0].op = sz;
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
  applyEmbeddings(0, 4, dims, computeBlockSize(4, dims), alpha, 1, embeddings,
                  x, y, 0, computeBlockSize(4, dims));
  for (int i = 0; i < computeBlockSize(4, dims); ++i) {
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
    EXPECT_FLOAT_EQ(y[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(y[i].im, expectedResult.im);
  }
  free(embeddings);

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
  applyEmbeddings(0, 4, dims, computeBlockSize(4, dims), alpha, 2, embeddings,
                  x, y, 0, computeBlockSize(4, dims));
  /* Make sure that dim of first space is 2.  Otherwise the expected
   * results are not calculated correctly. */
  EXPECT_EQ(dims[0], 2);
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
    EXPECT_FLOAT_EQ(y[i].re, expectedResult.re);
    EXPECT_FLOAT_EQ(y[i].im, expectedResult.im);
  }
  free(embeddings);

  mboElemOpDestroy(&sz);
  mboElemOpDestroy(&sp);

  free(x);
  free(y);
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
