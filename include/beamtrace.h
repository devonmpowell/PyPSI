#ifndef _BEAMTRACE_
#define _BEAMTRACE_

#include "psi.h" 
#include "geometry.h" 
#include "refine.h" 
#include "grid.h" 
#include "mesh.h" 
#include "rtree.h" 



// stores a beam: 3 vertices with corresponding 4-positions
// and 4-velocities 
typedef struct {

	psi_4vec pos; // 4-position
	psi_4vec vel; // 4-tangent vector
	psi_real alpha; // affine parameter

} psi_beam;


#define PSI_CLIGHT 1
#define PSI_METRIC_MINKOWSKI 0
#define PSI_METRIC_FLRW 1
#define PSI_METRIC_KERR 2

// structure to hold voxelization info
// remains the same regardless of dimensionality
typedef struct {

	// type flags and parameters for analytic metrics
	psi_int type;
	psi_real params[4];

	// these fields are used only for perturbed FLRW!
	psi_int snapnum;
	psi_dvec n;
	psi_rvec box;
	psi_real* phi;
	psi_real* gradphi;
} psi_metric;





void psi_beamtrace(psi_grid* grid, psi_mesh* mesh, psi_int bstep, psi_metric* metric, psi_real* outar, psi_dvec ardim);



#endif // _BEAMTRACE_
