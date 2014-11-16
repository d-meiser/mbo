#include <Integrator.h>
#include <stdlib.h>

struct RK4_ctx {
	struct MboAmplitude *k1, *k2, *k3, *k4;
};

void rk4_create(struct Integrator *self, MboGlobInd dim);
void rk4_destroy(struct Integrator *self);
void rk4_takeStep(struct Integrator *self, const struct MboAmplitude *x,
		  struct MboAmplitude *y, RHS f, void *ctx);
void rk4_advanceBeyond(struct Integrator *self, double t,
		       const struct MboAmplitude *x, struct MboAmplitude *y,
		       RHS f, void *ctx);
void rk4_advanceTo(struct Integrator *self, double t,
		   const struct MboAmplitude *x, struct MboAmplitude *y, RHS f,
		   void *ctx);

void integratorCreate(struct Integrator *integrator, MboGlobInd dim)
{
	integrator->ops.create = &rk4_create;
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
	self->data = ctx;
}

void rk4_destroy(struct Integrator* self) {
	struct RK4_ctx *ctx = (struct RK4_ctx *)self->data;
	if (ctx) {
		free(ctx->k1);
		free(ctx->k2);
		free(ctx->k3);
		free(ctx->k4);
		free(self->data);
	}
	self->ops.create = 0;
	self->ops.destroy = 0;
	self->ops.takeStep = 0;
	self->ops.advanceBeyond = 0;
	self->ops.advanceTo = 0;
	self->data = 0;
}

void rk4_takeStep(struct Integrator *self, const struct MboAmplitude *x,
		  struct MboAmplitude *y, RHS f, void *ctx)
{
}

void rk4_advanceBeyond(struct Integrator *self, double t,
		       const struct MboAmplitude *x, struct MboAmplitude *y,
		       RHS f, void *ctx)
{
}

void rk4_advanceTo(struct Integrator *self, double t,
		   const struct MboAmplitude *x, struct MboAmplitude *y, RHS f,
		   void *ctx)
{
}
