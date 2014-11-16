#include <Integrator.h>
#include <stdlib.h>

struct RK4_ctx {
	struct MboAmplitude *k1, *k2, *k3, *k4, *work;
};

void rk4_create(struct Integrator *self, MboGlobInd dim);
void rk4_destroy(struct Integrator *self);
void rk4_takeStep(struct Integrator *self, struct MboAmplitude *x, RHS f,
		  void *ctx);
void rk4_advanceBeyond(struct Integrator *self, double t,
		       struct MboAmplitude *x, RHS f, void *ctx);
void rk4_advanceTo(struct Integrator *self, double t, struct MboAmplitude *x,
		   RHS f, void *ctx);

void integratorCreate(struct Integrator *integrator, MboGlobInd dim)
{
	integrator->ops.create = &rk4_create;
	integrator->t = 0;
	integrator->dt = 1.0e-3;
	integrator->dim = dim;
	integrator->data = 0;
	integrator->ops.create(integrator, dim);
}

void integratorDestroy(struct Integrator* integrator)
{
	if (integrator->ops.destroy) {
		integrator->ops.destroy(integrator);
	}
}

void integratorSetTime(struct Integrator *integrator, double t)
{
	integrator->t = t;
}

double integratorGetTime(struct Integrator *integrator)
{
	return integrator->t;
}

void integratorTimeStepHint(struct Integrator *integrator, double dt)
{
	integrator->dt = dt;
}

void integratorTakeStep(struct Integrator *integrator, struct MboAmplitude *x,
			RHS f, void *ctx)
{
	integrator->ops.takeStep(integrator, x, f, ctx);
}

void integratorAdvanceBeyond(struct Integrator *integrator, double t,
			     struct MboAmplitude *x, RHS f, void *ctx)
{
	integrator->ops.advanceBeyond(integrator, t, x, f, ctx);
}

void integratorAdvanceTo(struct Integrator *integrator, double t,
			 struct MboAmplitude *x, RHS f, void *ctx)
{
	integrator->ops.advanceTo(integrator, t, x, f, ctx);
}

/* Implementation of RK4 integrator */

void rk4_create(struct Integrator *self, MboGlobInd dim)
{
	self->ops.destroy = &rk4_destroy;
	self->ops.takeStep = &rk4_takeStep;
	self->ops.advanceBeyond = &rk4_advanceBeyond;
	self->ops.advanceTo = &rk4_advanceTo;
	struct RK4_ctx *ctx = malloc(sizeof(*ctx));
	ctx->k1 = malloc(dim * sizeof(*ctx->k1));
	ctx->k2 = malloc(dim * sizeof(*ctx->k2));
	ctx->k3 = malloc(dim * sizeof(*ctx->k3));
	ctx->k4 = malloc(dim * sizeof(*ctx->k4));
	ctx->work = malloc(dim * sizeof(*ctx->work));
	self->data = ctx;
}

void rk4_destroy(struct Integrator* self) {
	struct RK4_ctx *ctx = (struct RK4_ctx *)self->data;
	if (ctx) {
		free(ctx->k1);
		free(ctx->k2);
		free(ctx->k3);
		free(ctx->k4);
		free(ctx->work);
		free(self->data);
	}
	self->ops.create = 0;
	self->ops.destroy = 0;
	self->ops.takeStep = 0;
	self->ops.advanceBeyond = 0;
	self->ops.advanceTo = 0;
	self->data = 0;
}

static void zaxpy(struct MboAmplitude *w, double alpha,
		  const struct MboAmplitude *x, const struct MboAmplitude *y,
		  MboGlobInd dim)
{
	MboGlobInd i;
	for (i = 0; i < dim; ++i) {
		w[i].re = alpha * x[i].re + y[i].re;
		w[i].im = alpha * x[i].im + y[i].im;
	}
}

void rk4_takeStep(struct Integrator *self, struct MboAmplitude *x, RHS f,
		  void *ctx)
{
	double prefactor;
	MboGlobInd i;

	struct RK4_ctx *rk4ctx = (struct RK4_ctx *)self->data;
	f(self->t, x, rk4ctx->k1, ctx);
	zaxpy(rk4ctx->work, 0.5 * self->dt, rk4ctx->k1, x, self->dim);
	f(self->t + 0.5 * self->dt, rk4ctx->work, rk4ctx->k2, ctx);
	zaxpy(rk4ctx->work, 0.5 * self->dt, rk4ctx->k2, x, self->dim);
	f(self->t + 0.5 * self->dt, rk4ctx->work, rk4ctx->k3, ctx);
	zaxpy(rk4ctx->work, self->dt, rk4ctx->k3, x, self->dim);
	f(self->t + self->dt, rk4ctx->work, rk4ctx->k4, ctx);
	prefactor = self->dt / 6.0;
	for (i = 0; i < self->dim; ++i) {
		x[i].re +=
		    prefactor * (rk4ctx->k1[i].re +
				 2.0 * (rk4ctx->k2[i].re + rk4ctx->k3[i].re) +
				 rk4ctx->k4[i].re);
		x[i].im +=
		    prefactor * (rk4ctx->k1[i].im +
				 2.0 * (rk4ctx->k2[i].im + rk4ctx->k3[i].im) +
				 rk4ctx->k4[i].im);
	}
	self->t += self->dt;
}

void rk4_advanceBeyond(struct Integrator *self, double t,
		       struct MboAmplitude *x, 
		       RHS f, void *ctx)
{
	while (self->t < t) {
		rk4_takeStep(self, x, f, ctx);
	}
}

void rk4_advanceTo(struct Integrator *self, double t,
		   struct MboAmplitude *x, RHS f,
		   void *ctx)
{
	double saveDt;
	while (self->t + self->dt < t) {
		rk4_takeStep(self, x, f, ctx);
	}
	saveDt = self->dt;
	self->dt = t - self->t;
	rk4_takeStep(self, x, f, ctx);
	self->dt = saveDt;
}
