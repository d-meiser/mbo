#include <SimpleTOpTraversal.h>

struct Tile tileIntersection(struct Tile t1, struct Tile t2)
{
	struct Tile intersection;
	intersection.rmin = (t1.rmin > t2.rmin) ? t1.rmin : t2.rmin;
	intersection.rmax = (t1.rmax < t2.rmax) ? t1.rmax : t2.rmax;
	intersection.cmin = (t1.cmin > t2.cmin) ? t1.cmin : t2.cmin;
	intersection.cmax = (t1.cmax < t2.cmax) ? t1.cmax : t2.cmax;
	return intersection;
}

int tileIsEmpty(const struct Tile *tile)
{
	return tile->rmin >= tile->rmax || tile->cmin >= tile->cmax;
}
