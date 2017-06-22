#ifndef _MESH_H_
#define _MESH_H_

#include "psi.h"

// structure to hold voxelization info
// remains the same regardless of dimensionality
typedef struct {
	char* type;
	psi_rvec *pos, *vel;
	psi_real *mass;
	psi_int npart;
	psi_int* connectivity;
} psi_mesh;


void psi_mesh_destroy(psi_mesh* mesh);


// mesh loaders (more to be added)


int peek_gadget2(psi_mesh* mesh, const char* filename);
int load_gadget2(psi_mesh* mesh, const char* filename);




#endif // _MESH_H_
