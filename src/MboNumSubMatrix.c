#include <MboNumSubMatrix.h>
#include <Tile.h>
#include <stdlib.h>

struct MboNumSubMatrix_t {
	MboNumOp op;
	struct Tile tile;
};

MboNumSubMatrix mboNumSubMatrixCreate(MboNumOp op, MboGlobInd rmin,
				      MboGlobInd rmax, MboGlobInd cmin,
				      MboGlobInd cmax)
{
	MboNumSubMatrix submat = malloc(sizeof(*submat));
	submat->op = op;
	submat->tile.rmin = rmin;
	submat->tile.rmax = rmax;
	submat->tile.cmin = cmin;
	submat->tile.cmax = cmax;
	return submat;
}

void mboNumSubMatrixDestroy(MboNumSubMatrix *m) 
{
	free(*m);
	*m = 0;
}

void mboNumSubMatrixSetTile(MboNumSubMatrix m, MboGlobInd rmin, MboGlobInd rmax,
			    MboGlobInd cmin, MboGlobInd cmax)
{
	m->tile.rmin = rmin;
	m->tile.rmax = rmax;
	m->tile.cmin = cmin;
	m->tile.cmax = cmax;
}
