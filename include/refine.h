/*
 *
 * refine.h 
 *
 * Copyright (c) 2013-2023 Devon M. Powell
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#ifndef _REFINE_H_ 
#define _REFINE_H_ 

#include "psi.h"
#include "mesh.h"

// a struct to wrap the refinement buffer
typedef struct {
	psi_rvec* pos;
	psi_rvec* vel;
	psi_real* mass;
	psi_int max_lvl;
	psi_real ptol2; // vtol2;
	//psi_int use_vel;
	psi_int num;
	psi_int capacity;
} psi_tet_buffer;

void psi_tet_buffer_init(psi_tet_buffer* buffer, psi_real tol, psi_int max_lvl);

void psi_tet_buffer_destroy(psi_tet_buffer* buffer);

void psi_tet_buffer_grow(psi_tet_buffer* buffer);

void psi_tet_buffer_refine(psi_tet_buffer* buffer, psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_int order);

#endif // _REFINE_H_ 
