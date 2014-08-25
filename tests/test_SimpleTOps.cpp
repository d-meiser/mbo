#include <gtest/gtest.h>
#include <SimpleTOp.h>

TEST(SimpleTOps, Kron) {
  struct SimpleTOp a, b, c;

  a.numFactors = 0;
  a.embeddings = 0;
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 0, &b, &c);
  EXPECT_EQ(c.numFactors, 0);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 0;
  a.embeddings = 0;
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 1;
  c.embeddings = (struct Embedding*)malloc(sizeof(*c.embeddings));
  c.embeddings[0].i = 0;
  mboElemOpCreate(&c.embeddings[0].op);
  kronSimpleTOps(&a, 0, &b, &c);
  EXPECT_EQ(c.numFactors, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 1;
  a.embeddings = (struct Embedding*)malloc(sizeof(*a.embeddings));
  a.embeddings[0].i = 0;
  mboElemOpCreate(&a.embeddings[0].op);
  b.numFactors = 0;
  b.embeddings = 0;
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 1, &b, &c);
  EXPECT_EQ(c.numFactors, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);

  a.numFactors = 1;
  a.embeddings = (struct Embedding*)malloc(sizeof(*a.embeddings));
  a.embeddings[0].i = 0;
  mboElemOpCreate(&a.embeddings[0].op);
  b.numFactors = 1;
  b.embeddings = (struct Embedding*)malloc(sizeof(*b.embeddings));
  b.embeddings[0].i = 0;
  mboElemOpCreate(&b.embeddings[0].op);
  c.numFactors = 0;
  c.embeddings = 0;
  kronSimpleTOps(&a, 1, &b, &c);
  EXPECT_EQ(c.numFactors, 2);
  EXPECT_EQ(c.embeddings[0].i, 0);
  EXPECT_EQ(c.embeddings[1].i, 1);
  destroySimpleTOp(&a);
  destroySimpleTOp(&b);
  destroySimpleTOp(&c);
}

