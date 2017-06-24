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
