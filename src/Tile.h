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
#ifndef TILE_H
#define TILE_H

#include <MboIndices.h>

#ifdef __cplusplus
extern "C" {
#endif

struct Tile
{
	MboGlobInd rmin;
	MboGlobInd rmax;
	MboGlobInd cmin;
	MboGlobInd cmax;
};

struct Tile tileIntersection(struct Tile t1, struct Tile t2);
int tileIsEmpty(const struct Tile *tile);
struct Tile tileDivide(struct Tile t1, struct Tile t2);
MboGlobInd numTilesContained(struct Tile t1, struct Tile t2);
void tileAdvance(MboGlobInd n, struct Tile *t);

#ifdef __cplusplus
}
#endif
#endif

