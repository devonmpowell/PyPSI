/*
 *
 *		psi.h
 *		
 *		Tools for conservative voxelization of finite elements.
 *		
 *		Devon Powell
 *		11 Feb 2016
 *
 */		

#ifndef _PSI_H_
#define _PSI_H_

// use python functions if compiled as a Python module
#ifdef PYMODULE
#include <Python.h>
#define psi_printf PySys_WriteStdout 
#define psi_malloc PyMem_Malloc
#define psi_realloc PyMem_Realloc 
#define psi_free PyMem_Free
// otherwise, just use the standard ones
#else
#define psi_printf printf
#define psi_malloc malloc
#define psi_realloc realloc
#define psi_free free
#endif

// standard includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/// TODO: this is hacky !!
// Make dimensionality dynamic!
#define PSI_NDIM 3



// primitive types 
// with floating-point precision as an option
typedef int16_t psi_short;
typedef int32_t psi_int;
typedef int64_t psi_long;
typedef double psi_real;

// union types for vectors and indexing
// along with other dimensionality-dependent things
typedef union {
	struct { psi_real x, y, z; };
	psi_real xyz[3];
} psi_rvec;
typedef union {
	struct { psi_int i, j, k; };
	psi_int ijk[3];
} psi_dvec;

// all other psi includes here, after primitive types 
//#include "grid.h"
//#include "geometry.h"
//#include "mesh.h"
//#include "skymap.h"

// sampling options
#define PSI_SAMPLING_VOLUME 0
#define PSI_SAMPLING_POINT 1
#define PSI_SAMPLING_CIC 2
#define PSI_SAMPLING_NGP 2

//void psi_voxels(psi_grid* grid, psi_mesh* mesh);


// user-exposed functions taking plain arrays 
//void psi_element_mesh(psi_rvec* pos, psi_rvec* vel, psi_real* mass, psi_int nelem, 
		//psi_int order, psi_real fractol, psi_dest_grid* grid);
//void psi_particle_mesh(psi_rvec* pos, psi_rvec* vel, psi_real* mass, psi_int nelem, psi_dest_grid* grid);
//void psi_element_block_from_grid(psi_rvec* pos_in, psi_rvec* vel_in, psi_real mtot, 
		//psi_dvec* ngrid, psi_int order, psi_dvec* start, psi_dvec* num, psi_rvec* pos_out, psi_rvec* vel_out, psi_real* mass_out);

//// user-exposed functions requiring an rtree as an argument
//#include "rtree.h"
//psi_int psi_sample_vdf(psi_rtree* rtree, psi_rvec samppos, psi_real reftol, psi_real* rhoout, psi_rvec* velout, psi_int capacity); 
//void psi_annihilate(psi_rtree* rtree, psi_real reftol, psi_dest_grid* grid); 


#endif // _PSI_H_
