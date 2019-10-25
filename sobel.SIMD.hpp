/*
  Copyright (C) 2017 Hoyoung Lee

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
#include <iostream>
#include <cstring>
#include <cstdint>

#include <cpixmap.hpp>

#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
# define MAX_VECTOR_SIZE 512
# include <vectorclass/vectorclass.h>
# if INSTRSET < 2
#  error "Unsupported x86-SIMD! Please comment USE_SIMD on!"
# endif
#elif defined(__GNUC__) && defined (__ARM_NEON__)
# include <arm_neon.h>
#else
# error "Undefined SIMD!"
#endif

inline void edgeHSobelKernel(cpixmap<uint8_t>& gray, cpixmap<int8_t>& dx)
{
  assert(gray.isMatched(dx));

  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint8_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int8_t *dxLine = dx.getLine(y, z);
      uint8_t *prevLine = gray3x3.getPrevLine();
      uint8_t *currLine = gray3x3.getCurrLine();
      uint8_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 32) {
	Vec32uc nwVec, /*nnVec,*/ neVec;
	nwVec.load(&prevLine[(int)x-1]), /*nnVec.load(&prevLine[(int)x+0]),*/ neVec.load(&prevLine[(int)x+1]);
	Vec32uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+0]);
	Vec32uc swVec, /*ssVec,*/ seVec;
	swVec.load(&nextLine[(int)x-1]), /*ssVec.load(&nextLine[(int)x+0]),*/ seVec.load(&nextLine[(int)x+1]);
	Vec32c dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16uc nwVec, /*nnVec,*/ neVec;
	nwVec.load(&prevLine[(int)x-1]), /*nnVec.load(&prevLine[(int)x+0]),*/ neVec.load(&prevLine[(int)x+1]);
	Vec16uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+0]);
	Vec16uc swVec, /*ssVec,*/ seVec;
	swVec.load(&nextLine[(int)x-1]), /*ssVec.load(&nextLine[(int)x+0]),*/ seVec.load(&nextLine[(int)x+1]);
	Vec16c dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	uint8x16_t nwVec, /*nnVec,*/ neVec;
	nwVec = vld1q_u8((const uint8_t *)&prevLine[(int)x-1]);
	//nnVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+1]);
	uint8x16_t wwVec, eeVec;
	wwVec = vld1q_u8((const uint8_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u8((const uint8_t *)&currLine[(int)x+1]);
	uint8x16_t swVec, /*ssVec,*/ seVec;
	swVec = vld1q_u8((const uint8_t *)&nextLine[(int)x-1]);
	//ssVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+1]);
	
	int8x16_t lsumVec, rsumVec;
	uint8x16_t tempVec;
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, nwVec, 3);
	tempVec = vsra_n_u8(tempVec, wwVec, 2);
	tempVec = vsra_n_u8(tempVec, swVec, 3);
	lsumVec = vreinterpretq_s8_u8(tempVec);
	
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, neVec, 3);
	tempVec = vsra_n_u8(tempVec, eeVec, 2);
	tempVec = vsra_n_u8(tempVec, seVec, 3);
	rsumVec = vreinterpretq_s8_u8(tempVec);
	
	int8x16_t dxVec = vsubq_s8(rsum, lsum);
	vst1q_s8((int8_t *)&dxLine[x], dxVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}

inline void edgeVSobelKernel(cpixmap<uint8_t>& gray, cpixmap<int8_t>& dy)
{
  assert(gray.isMatched(dy));

  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint8_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int8_t *dyLine = dy.getLine(y, z);
      uint8_t *prevLine = gray3x3.getPrevLine();
      //uint8_t *currLine = gray3x3.getCurrLine();
      uint8_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 32) {
	Vec32uc nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	/*
	Vec32uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+0]);
	*/
	Vec32uc swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec32c dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16uc nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	/*
	Vec16uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	*/
	Vec16uc swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec16c dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	uint8x16_t nwVec, nnVec, neVec;
	nwVec = vld1q_u8((const uint8_t *)&prevLine[(int)x-1]);
	nnVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+1]);
	/*
	uint8x16_t wwVec, eeVec;
	wwVec = vld1q_u8((const uint8_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u8((const uint8_t *)&currLine[(int)x+]);
	*/
	uint8x16_t swVec, ssVec, seVec;
	swVec = vld1q_u8((const uint8_t *)&nextLine[(int)x-1]);
	ssVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+1]);
	
	int8x16_t tsumVec, bsumVec;
	uint8x16_t tempVec;
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, nwVec, 3);
	tempVec = vsra_n_u8(tempVec, nnVec, 2);
	tempVec = vsra_n_u8(tempVec, neVec, 3);
	tsumVec = vreinterpretq_s8_u8(tempVec);
	
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, swVec, 3);
	tempVec = vsra_n_u8(tempVec, ssVec, 2);
	tempVec = vsra_n_u8(tempVec, seVec, 3);
	bsumVec = vreinterpretq_s8_u8(tempVec);
	
	int8x16_t dyVec = vsubq_s8(bsum, tsum);
	vst1q_s8((int8_t *)&dyLine[x], dyVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}

inline void edgeSobelKernel(cpixmap<uint8_t>& gray, cpixmap<int8_t>& dx, cpixmap<int8_t>& dy)
{
  assert(gray.isMatched(dx));
  assert(gray.isMatched(dy));
  
  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint8_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int8_t *dxLine = dx.getLine(y, z);
      int8_t *dyLine = dy.getLine(y, z);
      uint8_t *prevLine = gray3x3.getPrevLine();
      uint8_t *currLine = gray3x3.getCurrLine();
      uint8_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 32) {
	Vec32uc nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	Vec32uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+0]);
	Vec32uc swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec32c dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
	Vec32c dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16uc nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	Vec16uc wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	Vec16uc swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec16c dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
	Vec16c dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	uint8x16_t nwVec, nnVec, neVec;
	nwVec = vld1q_u8((const uint8_t *)&prevLine[(int)x-1]);
	nnVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u8((const uint8_t *)&prevLine[(int)x+1]);
	uint8x16_t wwVec, eeVec;
	wwVec = vld1q_u8((const uint8_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u8((const uint8_t *)&currLine[(int)x+]);
	uint8x16_t swVec, ssVec, seVec;
	swVec = vld1q_u8((const uint8_t *)&nextLine[(int)x-1]);
	ssVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u8((const uint8_t *)&nextLine[(int)x+1]);

	int8x16_t lsumVec, rsumVec;
	uint8x16_t tempVec;
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, nwVec, 3);
	tempVec = vsra_n_u8(tempVec, wwVec, 2);
	tempVec = vsra_n_u8(tempVec, swVec, 3);
	lsumVec = vreinterpretq_s8_u8(tempVec);
	
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, neVec, 3);
	tempVec = vsra_n_u8(tempVec, eeVec, 2);
	tempVec = vsra_n_u8(tempVec, seVec, 3);
	rsumVec = vreinterpretq_s8_u8(tempVec);
	
	int8x16_t dxVec = vsubq_s8(rsum, lsum);
	vst1q_s8((int8_t *)&dxLine[x], dxVec);

	int8x16_t tsumVec, bsumVec;
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, nwVec, 3);
	tempVec = vsra_n_u8(tempVec, nnVec, 2);
	tempVec = vsra_n_u8(tempVec, neVec, 3);
	tsumVec = vreinterpretq_s8_u8(tempVec);
	
	tempVec = vdupq_n_u8((uint8_t)0);
	tempVec = vsra_n_u8(tempVec, swVec, 3);
	tempVec = vsra_n_u8(tempVec, ssVec, 2);
	tempVec = vsra_n_u8(tempVec, seVec, 3);
	bsumVec = vreinterpretq_s8_u8(tempVec);
	
	int8x16_t dyVec = vsubq_s8(bsum, tsum);
	vst1q_s8((int8_t *)&dyLine[x], dyVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}

inline void edgeHSobelKernel(cpixmap<uint16_t>& gray, cpixmap<int16_t>& dx)
{
  assert(gray.isMatched(dx));

  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint16_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int16_t *dxLine = dx.getLine(y, z);
      uint16_t *prevLine = gray3x3.getPrevLine();
      uint16_t *currLine = gray3x3.getCurrLine();
      uint16_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16us nwVec, /*nnVec,*/ neVec;
	nwVec.load(&prevLine[(int)x-1]), /*nnVec.load(&prevLine[(int)x+0]),*/ neVec.load(&prevLine[(int)x+1]);
	Vec16us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	Vec16us swVec, /*ssVec,*/ seVec;
	swVec.load(&nextLine[(int)x-1]), /*ssVec.load(&nextLine[(int)x+0]),*/ seVec.load(&nextLine[(int)x+1]);
	Vec16s dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	Vec8us nwVec, /*nnVec,*/ neVec;
	nwVec.load(&prevLine[(int)x-1]), /*nnVec.load(&prevLine[(int)x+0]),*/ neVec.load(&prevLine[(int)x+1]);
	Vec8us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	Vec8us swVec, /*ssVec,*/ seVec;
	swVec.load(&nextLine[(int)x-1]), /*ssVec.load(&nextLine[(int)x+0]),*/ seVec.load(&nextLine[(int)x+1]);
	Vec8s dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	uint16x8_t nwVec, /*nnVec,*/ neVec;
	nwVec = vld1q_u16((const uint16_t *)&prevLine[(int)x-1]);
	//nnVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+1]);
	uint16x8_t wwVec, eeVec;
	wwVec = vld1q_u16((const uint16_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u16((const uint16_t *)&currLine[(int)x+1]);
	uint16x8_t swVec, /*ssVec,*/ seVec;
	swVec = vld1q_u16((const uint16_t *)&nextLine[(int)x-1]);
	//ssVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+1]);
	
	int16x8_t lsumVec, rsumVec;
	uint16x8_t tempVec;
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, nwVec, 3);
	tempVec = vsra_n_u16(tempVec, wwVec, 2);
	tempVec = vsra_n_u16(tempVec, swVec, 3);
	lsumVec = vreinterpretq_s16_u16(tempVec);
	
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, neVec, 3);
	tempVec = vsra_n_u16(tempVec, eeVec, 2);
	tempVec = vsra_n_u16(tempVec, seVec, 3);
	rsumVec = vreinterpretq_s16_u16(tempVec);
	
	int16x8_t dxVec = vsubq_s16(rsum, lsum);
	vst1q_s16((int16_t *)&dxLine[x], dxVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}

inline void edgeVSobelKernel(cpixmap<uint16_t>& gray, cpixmap<int16_t>& dy)
{
  assert(gray.isMatched(dy));
  

  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint16_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int16_t *dyLine = dy.getLine(y, z);
      uint16_t *prevLine = gray3x3.getPrevLine();
      //uint16_t *currLine = gray3x3.getCurrLine();
      uint16_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16us nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	/*
	Vec16us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	*/
	Vec16us swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec16s dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	Vec8us nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	/*
	Vec8us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	*/
	Vec8us swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec8s dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	uint16x8_t nwVec, nnVec, neVec;
	nwVec = vld1q_u16((const uint16_t *)&prevLine[(int)x-1]);
	nnVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+1]);
	/*
	uint16x8_t wwVec, eeVec;
	wwVec = vld1q_u16((const uint16_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u16((const uint16_t *)&currLine[(int)x+1]);
	*/
	uint16x8_t swVec, ssVec, seVec;
	swVec = vld1q_u16((const uint16_t *)&nextLine[(int)x-1]);
	ssVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+1]);
	
	int16x8_t tsumVec, bsumVec;
	uint16x8_t tempVec;
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, nwVec, 3);
	tempVec = vsra_n_u16(tempVec, nnVec, 2);
	tempVec = vsra_n_u16(tempVec, neVec, 3);
	tsumVec = vreinterpretq_s16_u16(tempVec);
	
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, swVec, 3);
	tempVec = vsra_n_u16(tempVec, ssVec, 2);
	tempVec = vsra_n_u16(tempVec, seVec, 3);
	bsumVec = vreinterpretq_s16_u16(tempVec);
	
	int16x8_t dyVec = vsubq_s16(bsum, tsum);
	vst1q_s16((int16_t *)&dyLine[x], dyVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}

inline void edgeSobelKernel(cpixmap<uint16_t>& gray, cpixmap<int16_t>& dx, cpixmap<int16_t>& dy)
{
  assert(gray.isMatched(dx));
  assert(gray.isMatched(dy));

  for (size_t z  = 0; z < gray.getBands(); ++z) {
    window3x3_frame<uint16_t> gray3x3(gray);
    gray3x3.draftFrame(gray, z);

    for (size_t y = 0; y < gray.getHeight(); ++y) {
      int16_t *dxLine = dx.getLine(y, z);
      int16_t *dyLine = dy.getLine(y, z);
      uint16_t *prevLine = gray3x3.getPrevLine();
      uint16_t *currLine = gray3x3.getCurrLine();
      uint16_t *nextLine = gray3x3.getNextLine();

      /*
      nwVec|nnVec|neVec
      -----+-----+-----
      wwVec|ooVec|eeVec
      -----+-----+-----
      swVec|ssVec|seVec
      */
#pragma omp parallel for
#if defined(__x86_64__) || defined(__i386__)
# if INSTRSET >= 8 // AVXx - 256bits
      for (size_t x = 0; x < gray.getWidth(); x += 16) {
	Vec16us nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	Vec16us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	Vec16us swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec16s dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
	Vec16s dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# elif INSTRSET >= 2 // SSE2 - 128bits
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	Vec8us nwVec, nnVec, neVec;
	nwVec.load(&prevLine[(int)x-1]), nnVec.load(&prevLine[(int)x+0]), neVec.load(&prevLine[(int)x+1]);
	Vec8us wwVec, eeVec;
	wwVec.load(&currLine[(int)x-1]), eeVec.load(&currLine[(int)x+1]);
	Vec8us swVec, ssVec, seVec;
	swVec.load(&nextLine[(int)x-1]), ssVec.load(&nextLine[(int)x+0]), seVec.load(&nextLine[(int)x+1]);
	Vec8s dxVec =
	  -(nwVec>>3) + (neVec>>3)
	  -(wwVec>>2) + (eeVec>>2)
	  -(swVec>>3) + (seVec>>3);
	dxVec.store(&dxLine[x]);
	Vec8s dyVec =
	  -(nwVec>>3) - (nnVec>>2) - (neVec>>3)
	  +(swVec>>3) + (ssVec>>2) + (seVec>>3);
	dyVec.store(&dyLine[x]);
      }
# endif
#elif defined(__ARM_NEON__)
      for (size_t x = 0; x < gray.getWidth(); x += 8) {
	uint16x8_t nwVec, nnVec, neVec;
	nwVec = vld1q_u16((const uint16_t *)&prevLine[(int)x-1]);
	nnVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+0]);
	neVec = vld1q_u16((const uint16_t *)&prevLine[(int)x+1]);
	uint16x8_t wwVec, eeVec;
	wwVec = vld1q_u16((const uint16_t *)&currLine[(int)x-1]);
	eeVec = vld1q_u16((const uint16_t *)&currLine[(int)x+1]);
	uint16x8_t swVec, ssVec, seVec;
	swVec = vld1q_u16((const uint16_t *)&nextLine[(int)x-1]);
	ssVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+0]);
	seVec = vld1q_u16((const uint16_t *)&nextLine[(int)x+1]);
	
	int16x8_t lsumVec, rsumVec;
	uint16x8_t tempVec;
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, nwVec, 3);
	tempVec = vsra_n_u16(tempVec, wwVec, 2);
	tempVec = vsra_n_u16(tempVec, swVec, 3);
	lsumVec = vreinterpretq_s16_u16(tempVec);
	
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, neVec, 3);
	tempVec = vsra_n_u16(tempVec, eeVec, 2);
	tempVec = vsra_n_u16(tempVec, seVec, 3);
	rsumVec = vreinterpretq_s16_u16(tempVec);
	
	int16x8_t dxVec = vsubq_s16(rsum, lsum);
	vst1q_s16((int16_t *)&dxLine[x], dxVec);

	int16x8_t tsumVec, bsumVec;
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, nwVec, 3);
	tempVec = vsra_n_u16(tempVec, nnVec, 2);
	tempVec = vsra_n_u16(tempVec, neVec, 3);
	tsumVec = vreinterpretq_s16_u16(tempVec);
	
	tempVec = vdupq_n_u16((uint16_t)0);
	tempVec = vsra_n_u16(tempVec, swVec, 3);
	tempVec = vsra_n_u16(tempVec, ssVec, 2);
	tempVec = vsra_n_u16(tempVec, seVec, 3);
	bsumVec = vreinterpretq_s16_u16(tempVec);
	
	int16x8_t dyVec = vsubq_s16(bsum, tsum);
	vst1q_s16((int16_t *)&dyLine[x], dyVec);
      }
#endif
      gray3x3.shiftFrame(gray, z);
    }
  }
}
