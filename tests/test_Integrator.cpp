#include <gtest/gtest.h>
#include <Integrator.h>

TEST(Integrator, Create) {
  struct Integrator integrator;
  integratorCreate(&integrator, 2);
  integratorDestroy(&integrator);
}

TEST(Integrator, SetTime) {
  struct Integrator integrator;
  integratorCreate(&integrator, 2);
  integratorSetTime(&integrator, 3.5);
  double t = integratorGetTime(&integrator);
  EXPECT_FLOAT_EQ(3.5, t);
  integratorDestroy(&integrator);
}

TEST(Integrator, TimeStepHint) {
  struct Integrator integrator;
  integratorCreate(&integrator, 2);
  integratorTimeStepHint(&integrator, 1.0e-6);
  integratorDestroy(&integrator);
}

