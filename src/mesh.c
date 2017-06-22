#include "mesh.h"

void psi_mesh_destroy(psi_mesh* mesh) {

	psi_free(mesh->pos);
	psi_free(mesh->vel);
	psi_free(mesh->mass);
	psi_free(mesh->connectivity);
}

typedef struct {
	int npart[6];
	double mass[6];
	double time;   // scale factor
	double redshift;                    
	int flag_sfr;                        
	int flag_feedback;                   
	unsigned int npartTotal[6];          
	int flag_cooling;                    
	int num_files;                       
	double BoxSize;                      
	double Omega0;                       
	double OmegaLambda;                  
	double HubbleParam;                  
	int flag_stellarage;                 
	int flag_metals;                     
	unsigned int npartTotalHighWord[6];
	int  flag_entropy_instead_u;         
	char fill[52];
	int max_level;
	int base_level;	
} gadget2header; 

int peek_gadget2(psi_mesh* mesh, const char* filename) {

	int blksize, e;
	gadget2header header;
	FILE* f = fopen(filename, "r");
	if(!f) {
		psi_printf("Failed to open %s.\n", filename);
		return 0;
	} 
	
	e = fread(&blksize, sizeof(int), 1, f);
	e = fread(&header, sizeof(gadget2header), 1, f);
	e = fread(&blksize, sizeof(int), 1, f);
	mesh->npart = header.npart[1]; 
	fclose(f);
	return 1;
}

int load_gadget2(psi_mesh* mesh, const char* filename) {

	int blksize, p, e;
	gadget2header header;
	FILE* f = fopen(filename, "r");
	if(!f) {
		psi_printf("Failed to open %s.\n", filename);
		return 0;
	} 
	e = fread(&blksize, sizeof(int), 1, f);
	e = fread(&header, sizeof(gadget2header), 1, f);
	e = fread(&blksize, sizeof(int), 1, f);

	// allocate mesh storage
	mesh->npart = header.npart[1]; 
	float* tpos = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	float* tvel = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	int* tid = (int*) psi_malloc(mesh->npart*sizeof(psi_int));
	
	// read in position data
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tpos, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==3*mesh->npart*sizeof(float));
	
	// read velocities
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tvel, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==3*mesh->npart*sizeof(float));
	
	// read ids
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==mesh->npart*sizeof(int));
	e = fread(tid, sizeof(int), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	assert(blksize==mesh->npart*sizeof(int));
	
	// close the input file
	fclose(f);
	
	// make arrays for vertices
	// NOTE: assumes the arrays have been malloced!
	for(p = 0; p < mesh->npart; ++p) {
		mesh->pos[tid[p]].x = tpos[3*p+0];
		mesh->pos[tid[p]].y = tpos[3*p+1];
		mesh->pos[tid[p]].z = tpos[3*p+2];
		mesh->vel[tid[p]].x = tvel[3*p+0];
		mesh->vel[tid[p]].y = tvel[3*p+1];
		mesh->vel[tid[p]].z = tvel[3*p+2];
	}

	psi_free(tpos);
	psi_free(tvel);
	psi_free(tid);
	//psi_printf("...done.\n");
	return 1;
}
