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
#include <Utilities.h>

MboGlobInd computeBlockSize(int N, MboLocInd *dims)
{
	int i;
	MboGlobInd blockSize = 1;
	for (i = 0; i < N; ++i) {
		blockSize *= (MboGlobInd)dims[i];
	}
	return blockSize;
}


