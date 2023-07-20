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
