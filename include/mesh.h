/*
 *
 * mesh.h 
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



#ifndef _MESH_H_
#define _MESH_H_

#include "psi.h"
#include "geometry.h"

#define PSI_MESH_SIMPLEX 4
#define PSI_MESH_LINEAR 8
#define PSI_MESH_QUADRATIC 27 


// structure to hold voxelization info
// remains the same regardless of dimensionality
typedef struct {
	psi_int elemtype;
	psi_int dim;
	psi_rvec *pos, *vel;
	psi_real *mass;
	psi_int npart;
	psi_int nelem;
	psi_int* connectivity;
	psi_rvec box[2];
	psi_int periodic;
} psi_mesh;

void psi_mesh_destroy(psi_mesh* mesh);

// mesh loaders (more to be added)

#define PSI_MESH_GADGET2 0
int peek_gadget2(psi_mesh* mesh, const char* filename);
int load_gadget2(psi_mesh* mesh, const char* filename);

#define PSI_MESH_GEVOLUTION 1
int load_gevolution(psi_mesh* mesh, const char* filename);

#endif // _MESH_H_
