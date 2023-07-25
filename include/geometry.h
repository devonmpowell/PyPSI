#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "psi.h"
#include "grid.h"

// numeric constants
#define PI 3.14159265359
#define TWO_PI 6.28318530718 
#define FOUR_PI 12.5663706144 
#define ONE_THIRD 0.3333333333333333333333333333333333
#define ONE_SIXTH 0.1666666666666666666666666666666666
#define ONE_TWLTH 0.0833333333333333333333333333333333


typedef struct {
	psi_rvec n;
	psi_real d;
} psi_plane;

// structure to hold a polyhedron,
// including mass coordinates
#define POLYSZ 64 // TODO: fix this number better. Should be 16, accd to Euler Characteristic with 10 faces
typedef struct {
	psi_rvec pos, q;
	psi_short pnbrs[PSI_NDIM];
} psi_vertex;
typedef struct {
	psi_int nverts;
	psi_dvec ibox[2];
	psi_vertex verts[POLYSZ];
} psi_poly; 



#define psi_dot3(va, vb) (va.x*vb.x + va.y*vb.y + va.z*vb.z)
#define psi_cross3(v,u,w)          /* CROSS Vector Product */              \
{                                                                       \
    (v).x = (u).y*(w).z - (u).z*(w).y;                             \
    (v).y = (u).z*(w).x - (u).x*(w).z;                             \
    (v).z = (u).x*(w).y - (u).y*(w).x;                             \
}

psi_real psi_omega3(psi_rvec v1, psi_rvec v2, psi_rvec v3);

psi_real psi_barycentric(psi_rvec* pos, psi_rvec samppos, psi_real* bcoords);

void psi_clip_reduce_tet(psi_rvec* pos, psi_plane* clip_planes, psi_int nclip, psi_real* moments);

void psi_voxelize_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_grid* grid);

void psi_voxelize_annihilation(psi_rvec* pos0, psi_rvec* vel0, psi_real mass0, 
		psi_rvec* pos1, psi_rvec* vel1, psi_real mass1, psi_rvec* mbox, psi_grid* grid);

//void psi_point_sample_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_dest_grid* dest_grid);


#endif // _GEOMETRY_H_
