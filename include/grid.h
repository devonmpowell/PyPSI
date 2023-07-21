/*
 *
 * grid.h 
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


#ifndef _GRID_H_
#define _GRID_H_

#include "psi.h"

#define PSI_GRID_CART 0
#define PSI_GRID_HPRING 1

// field tags
#define PSI_GRID_M 0 // weight (usually mass)
#define PSI_GRID_X 1 // position
#define PSI_GRID_V 2 // velocity
#define PSI_GRID_XX 3
#define PSI_GRID_XV 4
#define PSI_GRID_VV 5

// structure to hold voxelization info
// remains the same regardless of dimensionality
typedef struct {
	psi_int type;
	psi_int dim; 
	psi_real* fields[16];
	psi_rvec window[2];
	psi_dvec n;
	psi_rvec d;
} psi_grid;


// functions
//psi_int psi_grid_get_cell_geometry(psi_grid* grid, psi_dvec grind, psi_int bstep, psi_rvec* boundary, psi_rvec* centeri, psi_real* vol);

#endif // _GRID_H_
