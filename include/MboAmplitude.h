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
/**
 * @file MboAmplitude.h
 * @brief Definition of the complex amplitude type
 * */
#ifndef MBO_AMPLITUDE_H
#define MBO_AMPLITUDE_H

#ifdef __cplusplus
extern "C" {
#endif

struct MboAmplitude {
	double re;
	double im;
};

#ifdef __cplusplus
}
#endif
#endif
