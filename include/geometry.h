#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "psi.h"

// numeric constants
#define PI 3.14159265359
#define FOUR_PI 12.5663706144 
#define ONE_THIRD 0.3333333333333333333333333333333333
#define ONE_SIXTH 0.1666666666666666666666666666666666
#define ONE_TWLTH 0.0833333333333333333333333333333333

psi_real psi_omega3(psi_rvec v1, psi_rvec v2, psi_rvec v3);

//void psi_voxelize_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_dest_grid* dest_grid);

//void psi_point_sample_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_dest_grid* dest_grid);

//psi_real psi_barycentric(psi_rvec* pos, psi_rvec samppos, psi_real* bcoords);

//void psi_voxelize_annihilation(psi_rvec* pos0, psi_rvec* vel0, psi_real mass0, psi_rvec* rbox0, psi_rvec* pos1, psi_rvec* vel1, psi_real mass1, psi_rvec* rbox1, psi_dest_grid* dest_grid);


#endif // _GEOMETRY_H_
