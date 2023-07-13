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
