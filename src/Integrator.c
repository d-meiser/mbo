#include <Integrator.h>

void integratorCreate(struct Integrator* integrator)
{
	integrator->ops.create = 0;
	integrator->ops.destroy= 0;
	integrator->ops.takeStep = 0;
	integrator->ops.advanceBeyond= 0;
	integrator->ops.advanceTo = 0;
	integrator->data = 0;
}

void integratorDestroy(struct Integrator* integrator)
{
	if (integrator->ops.destroy) {
		integrator->ops.destroy(integrator);
	}
}
