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

#ifdef __cplusplus
extern "C" {
#endif

struct Integrator;

struct IntegratorOps {
	void (*create)(struct Integrator* self, MboGlobInd dim);
	void (*takeStep)(struct Integrator *self, const double *x,
			 double *y,
			 void (*f)(const double *x, double *y, void *ctx),
			 void *ctx);
	void (*advanceBeyond)(struct Integrator *self, double t,
			      const double *x, double *y,
			      void (*f)(const double *x, double *y, void *ctx),
			      void *ctx);
	void (*advanceTo)(struct Integrator *self, double t,
			  const double *x, double *y,
			  void (*f)(const double *x, double *y, void *ctx),
			  void *ctx);
	void (*destroy)(struct Integrator *self);
};

struct Integrator {
	struct IntegratorOps ops;
	double t;
	double dt;
	void *data;
};

void integratorCreate(struct Integrator* integrator);
void integratorDestroy(struct Integrator* integrator);

#ifdef __cplusplus
}
#endif

#endif

