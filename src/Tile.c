#include <Tile.h>


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

struct Tile tileDivide(struct Tile t1, struct Tile t2)
{
  struct Tile quotient;
  quotient.rmin = (t1.rmin - t2.rmin) / (t2.rmax - t2.rmin);
  quotient.rmax = (t1.rmax - t2.rmin) / (t2.rmax - t2.rmin);
  quotient.cmin = (t1.cmin - t2.cmin) / (t2.cmax - t2.cmin);
  quotient.cmax = (t1.cmax - t2.cmin) / (t2.cmax - t2.cmin);
  return quotient;
}

MboGlobInd numTilesContained(struct Tile t1, struct Tile t2)
{
  struct Tile quotient = tileDivide(t1, t2);
  MboGlobInd numRows = quotient.rmax - quotient.rmin;
  MboGlobInd numCols = quotient.cmax - quotient.cmin;
  return (numRows < numCols) ? numRows : numCols;
}

void tileAdvance(MboGlobInd n, struct Tile *t)
{
  MboGlobInd dr = t->rmax - t->rmin;
  MboGlobInd dc = t->cmax - t->cmin;
  t->rmin += n * dr;
  t->rmax += n * dr;
  t->cmin += n * dc;
  t->cmax += n * dc;
}
