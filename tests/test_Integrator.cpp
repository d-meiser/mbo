#include <gtest/gtest.h>
#include <Integrator.h>

TEST(Integrator, Create) {
  struct Integrator integrator;
  integratorCreate(&integrator);
  integratorDestroy(&integrator);
}
