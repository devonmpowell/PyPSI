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

// standard includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// use python functions if compiled as a Python module
#ifdef PYMODULE
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define psi_printf PySys_WriteStdout 
#define psi_malloc PyMem_Malloc
#define psi_realloc PyMem_Realloc 
#define psi_free PyMem_Free
#else
// otherwise, just use the standard ones
#define psi_printf printf
#define psi_malloc malloc
#define psi_realloc realloc
#define psi_free free
#endif

#define psi_assert(x) { \
	if(!(x)) { \
		psi_printf("Assert failed! %s, line %d.\n", __FILE__, __LINE__); \
	} \
}

/// TODO: this is hacky !!
// Make dimensionality dynamic!
#define PSI_NDIM 3
#define PRINT_EVERY 1024

// sampling options
#define PSI_SAMPLING_VOLUME 0
#define PSI_SAMPLING_POINT 1
#define PSI_SAMPLING_CIC 2
#define PSI_SAMPLING_NGP 2

// mode - decay or annihilation
#define PSI_MODE_DENSITY 0
#define PSI_MODE_ANNIHILATION 1


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
	struct { psi_real t, x, y, z; };
	psi_real txyz[4];
} psi_4vec;
typedef union {
	struct { psi_int i, j, k; };
	psi_int ijk[3];
} psi_dvec;


#endif // _PSI_H_
