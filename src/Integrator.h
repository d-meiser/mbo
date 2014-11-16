/*
Copyright 2014 Dominic Meiser

This file is part of mbo.

mbo is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

mbo is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with mbo.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <MboIndices.h>
#include <MboAmplitude.h>

#ifdef __cplusplus
extern "C" {
#endif

struct Integrator;
typedef void (*RHS)(double, const struct MboAmplitude *, struct MboAmplitude *,
		 void *);

struct IntegratorOps {
	void (*create)(struct Integrator *self, MboGlobInd dim);
	void (*takeStep)(struct Integrator *self, struct MboAmplitude *x, RHS f,
			 void *ctx);
	void (*advanceBeyond)(struct Integrator *self, double t,
			      struct MboAmplitude *x, RHS f, void *ctx);
	void (*advanceTo)(struct Integrator *self, double t,
			  struct MboAmplitude *x, RHS f, void *ctx);
	void (*destroy)(struct Integrator *self);
};

struct Integrator {
	struct IntegratorOps ops;
	double t;
	double dt;
	MboGlobInd dim;
	void *data;
};

void integratorCreate(struct Integrator* integrator, MboGlobInd dim);
void integratorDestroy(struct Integrator* integrator);
void integratorSetTime(struct Integrator* integrator, double t);
double integratorGetTime(struct Integrator* integrator);
void integratorTimeStepHint(struct Integrator* integrator, double dt);
void integratorTakeStep(struct Integrator *integrator, struct MboAmplitude *x,
			RHS f, void *ctx);
void integratorAdvanceBeyond(struct Integrator *integrator, double t,
			     struct MboAmplitude *x, RHS f, void *ctx);

#ifdef __cplusplus
}
#endif

#endif

