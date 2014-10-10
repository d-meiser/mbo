#include <stdlib.h>
#include <MboNumSubMatrix.h>
#include <MboNumOpPrivate.h>
#include <SimpleTOp.h>
#include <Tile.h>

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

MBO_STATUS mboNumSubMatrixMatVec(struct MboAmplitude alpha, MboNumSubMatrix m,
				 struct MboAmplitude *x,
				 struct MboAmplitude beta,
				 struct MboAmplitude *y)
{
	MboGlobInd r;
	struct MboAmplitude tmp;
	int i;
	MBO_STATUS err;

	for (r = 0; r < m->tile.rmax - m->tile.rmin; ++r) {
		tmp.re = beta.re * y[r].re - beta.im * y[r].im;
		tmp.im = beta.re * y[r].im + beta.im * y[r].re;
		y[r].re = tmp.re;
		y[r].im = tmp.im;
	}
	for (i = 0; i < m->op->numTerms; ++i) {
		err = applySimpleTOp(m->op->space, alpha, m->op->sum + i, x, y,
				     &m->tile);
		if (err != MBO_SUCCESS) return err;
	}
	return MBO_SUCCESS;
}
