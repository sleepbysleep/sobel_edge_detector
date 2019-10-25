/*
  Copyright (C) 2014 Hoyoung Lee

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <cassert>
#include <limits>
#include <iostream>
#include <cstring>
//#include <float.h>
#include <type_traits>

#include <cpixmap.hpp>
#include <cchunk.hpp>

/* Shift-operated Kernel alternative to normal kernel
  4 3 4                            1 2 1
  3 2 3 => FIR coefficients 1/16 * 2 ? 2
  4 3 4                            1 2 1
*/

#if !defined(USE_SIMD)

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, void>::type
edgeHSobelKernel(cpixmap<T>& gray, cpixmap<typename std::make_signed<T>::type>& dx)
{
  assert(gray.isMatched(dx));
  
  using signed_T = typename std::make_signed<T>::type;

  for (size_t z = 0; z < gray.getBands(); ++z) {
    //window3x3_frame<typename std::make_unsigned<T>::type> gray3x3(gray);
    window3x3_frame<T> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      signed_T *dxLine = dx.getLine(y, z);
#pragma omp parallel for
      for (size_t x = 0; x < gray.getWidth(); ++x) {
	dxLine[x] =
	  -(gray3x3(y-1, x-1)>>3) + (gray3x3(y-1, x+1)>>3)
	  -(gray3x3(  y, x-1)>>2) + (gray3x3(  y, x+1)>>2)
	  -(gray3x3(y+1, x-1)>>3) + (gray3x3(y+1, x+1)>>3);
      }
      gray3x3.shiftFrame(gray, z);
    }
  }
}

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, void>::type
edgeVSobelKernel(cpixmap<T>& gray, cpixmap<typename std::make_signed<T>::type>& dy)
{
  assert(gray.isMatched(dy));
  
  using signed_T = typename std::make_signed<T>::type;
  
  for (size_t z = 0; z < gray.getBands(); ++z) {
    window3x3_frame<T> gray3x3(gray);
    gray3x3.draftFrame(gray, z);
    for (size_t y = 0; y < gray.getHeight(); ++y) {
      signed_T *dyLine = dy.getLine(y, z);
#pragma omp parallel for
      for (size_t x = 0; x < gray.getWidth(); ++x) {
	dyLine[x] =
	  -(gray3x3(y-1, x-1)>>3) - (gray3x3(y-1, x)>>2) - (gray3x3(y-1, x+1)>>3)
	  +(gray3x3(y+1, x-1)>>3) + (gray3x3(y+1, x)>>2) + (gray3x3(y+1, x+1)>>3);
      }
      gray3x3.shiftFrame(gray, z);
    }
  }
}

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, void>::type
edgeSobelKernel(cpixmap<T>& gray,
		cpixmap<typename std::make_signed<T>::type>& dx,
		cpixmap<typename std::make_signed<T>::type>& dy)
{
  assert(gray.isMatched(dx));
  assert(gray.isMatched(dy));

  using signed_T = typename std::make_signed<T>::type;
  
  for (size_t z = 0; z < gray.getBands(); ++z) {
    window3x3_frame<T> gray3x3(gray);
    gray3x3.draftFrame(gray, z);
    for (size_t y = 0; y < gray.getHeight(); ++y) {
      signed_T *dxLine = dx.getLine(y, z);
      signed_T *dyLine = dy.getLine(y, z);
#pragma omp parallel for
      for (size_t x = 0; x < gray.getWidth(); ++x) {
	dxLine[x] =
	  -(gray3x3(y-1, x-1)>>3) + (gray3x3(y-1, x+1)>>3)
	  -(gray3x3(  y, x-1)>>2) + (gray3x3(  y, x+1)>>2)
	  -(gray3x3(y+1, x-1)>>3) + (gray3x3(y+1, x+1)>>3);
	dyLine[x] =
	  -(gray3x3(y-1, x-1)>>3) - (gray3x3(y-1, x)>>2) - (gray3x3(y-1, x+1)>>3)
	  +(gray3x3(y+1, x-1)>>3) + (gray3x3(y+1, x)>>2) + (gray3x3(y+1, x+1)>>3);
      }
      gray3x3.shiftFrame(gray, z);
    }
  }
}

#else
# include "sobel.SIMD.hpp"
#endif
