#ifndef _SKYMAP_
#define _SKYMAP_

#include "psi.h" 
#include "geometry.h" 
#include "refine.h" 
#include "grid.h" 
#include "mesh.h" 
#include "rtree.h" 

//#define PI_OVER_THREE 1.0471975512
//#define PI 3.14159265359 
//#define TWO_PI 6.28318530718 
//#define SQRT_TWO_PI 2.50662827463 
//#define FOUR_PI 12.5663706144
//#define ONE_OVER_FOUR_PI 0.07957747154 
//#define TO_RADIANS 0.01745329251
//#define NPOINT_SIGMA 2048 

void psi_skymap(psi_grid* grid, psi_mesh* mesh, psi_int bstep);


// stores a beam: 3 vertices with corresponding 4-positions
// and 4-velocities 
typedef struct {

	psi_real pos[3][4];
	psi_real vel[4][4];
	psi_real flux;

} psi_beam;



#endif // _SKYMAP_
