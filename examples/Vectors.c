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
/*
\par
Illustrates the use of vectors in MBO.
*/
#include <MboVec.h>
#include <MboAmplitude.h>

int main()
{
	struct MboAmplitude alpha = {2.0, 3.0};
	long d = 10l;
	int i;
	MboVec x;
	MboVec y;
	mboVecCreate(d, &x);
	mboVecCreate(d, &y);

	for (i = 0; i < 10; ++i) {
		mboVecAXPY(&alpha, x, y);
	}

	mboVecDestroy(&x);
	mboVecDestroy(&y);
	return 0;
}

