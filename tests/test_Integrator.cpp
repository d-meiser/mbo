#include <gtest/gtest.h>
#include <cmath>
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

struct DecayCtx {
  double gamma;
};
void exponentialDecay(double t, const struct MboAmplitude* x,
                      struct MboAmplitude* y, void* ctx) {
  struct DecayCtx* decCtx = (struct DecayCtx*)ctx;
  y[0].re = -decCtx->gamma * x[0].re;
  y[0].im = -decCtx->gamma * x[0].im;
}

TEST(Integrator, TakeStep) {
  struct Integrator integrator;
  integratorCreate(&integrator, 1);
  double dt = 1.0e-3;
  integratorTimeStepHint(&integrator, dt);
  struct MboAmplitude x;
  x.re = 1.0;
  x.im = 0.0;
  struct MboAmplitude y;
  struct DecayCtx ctx;
  ctx.gamma = 1.0;
  integratorTakeStep(&integrator, &x, &y, &exponentialDecay, &ctx);
  EXPECT_FLOAT_EQ(x.re * exp(-ctx.gamma * dt), y.re);
  EXPECT_FLOAT_EQ(0, y.im);
  integratorDestroy(&integrator);
}

