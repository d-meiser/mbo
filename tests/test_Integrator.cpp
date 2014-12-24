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
  EXPECT_DOUBLE_EQ(3.5, t);
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
  struct DecayCtx ctx;
  ctx.gamma = 1.0;
  integratorTakeStep(&integrator, &x, &exponentialDecay, &ctx);
  EXPECT_DOUBLE_EQ(exp(-ctx.gamma * dt), x.re);
  EXPECT_DOUBLE_EQ(0, x.im);
  integratorDestroy(&integrator);
}

TEST(Integrator, AdvanceBeyond) {
  struct Integrator integrator;
  integratorCreate(&integrator, 1);
  double dt = 1.0e-3;
  integratorTimeStepHint(&integrator, dt);
  struct MboAmplitude x;
  x.re = 1.0;
  x.im = 0.0;
  struct DecayCtx ctx;
  ctx.gamma = 1.0;
  integratorAdvanceBeyond(&integrator, 0.3, &x, &exponentialDecay, &ctx);
  double finalTime = integratorGetTime(&integrator);
  EXPECT_LE(0.3, finalTime);
  EXPECT_LE(std::abs(exp(-finalTime * ctx.gamma) - x.re), 1.0e-8);
  integratorDestroy(&integrator);
}

TEST(Integrator, AdvanceTo) {
  struct Integrator integrator;
  integratorCreate(&integrator, 1);
  double dt = 1.0e-3;
  integratorTimeStepHint(&integrator, dt);
  struct MboAmplitude x;
  x.re = 1.0;
  x.im = 0.0;
  struct DecayCtx ctx;
  ctx.gamma = 1.0;
  double targetTime = 0.33458;
  integratorAdvanceTo(&integrator, targetTime, &x, &exponentialDecay, &ctx);
  double finalTime = integratorGetTime(&integrator);
  EXPECT_DOUBLE_EQ(targetTime, finalTime);
  EXPECT_LE(std::abs(exp(-finalTime * ctx.gamma) - x.re), 1.0e-8);
  integratorDestroy(&integrator);
}

