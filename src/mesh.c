#include "mesh.h"





int psi_verts_per_elem(psi_int elemtype) {
	//const static psi_int vpere[3] = {4, 8, 27};
	return elemtype;
}



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

	int blksize, e, i;
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
	mesh->nelem = header.npart[1];
	mesh->periodic = 1;
	mesh->elemtype = PSI_MESH_LINEAR;
	mesh->dim = 3;
	for(i = 0; i < 3; ++i) {
		mesh->box[0].xyz[i] = 0.0;
		mesh->box[1].xyz[i] = header.BoxSize;
	}
	fclose(f);
	return 1;
}

int load_gadget2(psi_mesh* mesh, const char* filename) {

	int blksize, p, e, nside, ii, jj, kk, i, j, k;
	int elemind, locind, vertind;
	int load_mass;
	gadget2header header;
	psi_printf("Loading %s...\n", filename);
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
	mesh->nelem = header.npart[1];
	load_mass = (header.mass[1] == 0);
	mesh->periodic = 1;
	mesh->elemtype = PSI_MESH_LINEAR;
	mesh->dim = 3;
	for(i = 0; i < 3; ++i) {
		mesh->box[0].xyz[i] = 0.0;
		mesh->box[1].xyz[i] = header.BoxSize;
	}
	float* tpos = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	float* tvel = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	int* tid = (int*) psi_malloc(mesh->npart*sizeof(psi_int));

	printf("Hubble = %f\n", header.HubbleParam);
	printf("Box = %f\n", header.BoxSize);
	printf("load_mass = %d, mass = %f\n", load_mass, header.mass[1]);
	printf("n_part = %d\n", mesh->npart);

	// read in position data
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tpos, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));

	// read velocities
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tvel, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));

	// read ids
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==mesh->npart*sizeof(int));
	e = fread(tid, sizeof(int), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==mesh->npart*sizeof(int));

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
		mesh->mass[tid[p]] = header.mass[1]/mesh->npart;

	}

	psi_free(tpos);
	psi_free(tvel);
	psi_free(tid);

	// now build the mesh connectivity
	// trilinear elements naturally
	nside = floor(pow(mesh->npart+0.5, ONE_THIRD));
	for(i = 0; i < nside; ++i)
	for(j = 0; j < nside; ++j)
	for(k = 0; k < nside; ++k) {
		elemind = nside*nside*i + nside*j + k;
		for(ii = 0; ii < 2; ++ii)
		for(jj = 0; jj < 2; ++jj)
		for(kk = 0; kk < 2; ++kk) {
			locind = 4*ii + 2*jj + kk;
			vertind = nside*nside*((i+ii)%nside)
				+ nside*((j+jj)%nside) + ((k+kk)%nside);
			mesh->connectivity[8*elemind+locind] = vertind;
		}
	}
	psi_printf("...done.\n");
	return 1;
}


int load_gevolution(psi_mesh* mesh, const char* filename) {

	int blksize, p, e, nside, ii, jj, kk, i, j, k;
	int elemind, locind, vertind;
	int load_mass;
	gadget2header header;
	psi_printf("Loading %s...\n", filename);
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
	mesh->nelem = header.npart[1];
	load_mass = (header.mass[1] == 0);
	mesh->periodic = 1;
	mesh->elemtype = PSI_MESH_LINEAR;
	mesh->dim = 3;
	for(i = 0; i < 3; ++i) {
		mesh->box[0].xyz[i] = 0.0;
		mesh->box[1].xyz[i] = header.BoxSize;
	}
	float* tpos = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	float* tvel = (float*) psi_malloc(mesh->npart*sizeof(psi_rvec));
	long* tid = (long*) psi_malloc(mesh->npart*sizeof(psi_long));

	printf("Hubble = %f\n", header.HubbleParam);
	printf("Box = %f\n", header.BoxSize);
	printf("load_mass = %d, mass = %f\n", load_mass, header.mass[1]);
	printf("n_part = %d\n", mesh->npart);

	// read in position data
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tpos, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));

	// read velocities
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));
	e = fread(tvel, 3*sizeof(float), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==3*mesh->npart*sizeof(float));

	// read particle ids
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==mesh->npart*sizeof(long));
	e = fread(tid, sizeof(long), mesh->npart, f);
	e = fread(&blksize, sizeof(int), 1, f);
	psi_assert(blksize==mesh->npart*sizeof(long));

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
		mesh->mass[tid[p]] = header.mass[1]/mesh->npart;

	}

	psi_free(tpos);
	psi_free(tvel);
	psi_free(tid);

	// now build the mesh connectivity
	// trilinear elements naturally
	nside = floor(pow(mesh->npart+0.5, ONE_THIRD));
	for(i = 0; i < nside; ++i)
	for(j = 0; j < nside; ++j)
	for(k = 0; k < nside; ++k) {
		elemind = nside*nside*i + nside*j + k;
		for(ii = 0; ii < 2; ++ii)
		for(jj = 0; jj < 2; ++jj)
		for(kk = 0; kk < 2; ++kk) {
			locind = 4*ii + 2*jj + kk;
			vertind = nside*nside*((i+ii)%nside)
				+ nside*((j+jj)%nside) + ((k+kk)%nside);
			mesh->connectivity[8*elemind+locind] = vertind;
		}
	}
#if 1
	nside = floor(pow(mesh->npart+0.5, ONE_THIRD));

#define PPB 64
	for(p = 0; p < mesh->npart; ++p) {

		psi_int block = p/PPB;
		psi_int bind = p%PPB;
		psi_int bi = bind/16;
		psi_int bj = (bind - 16*bi)/4;
		psi_int bk = (bind - 16*bi - 4*bj);

		if(bi > 2 || bj > 2 || bk > 2)
			continue;

		for(ii = 0; ii < 2; ++ii)
		for(jj = 0; jj < 2; ++jj)
		for(kk = 0; kk < 2; ++kk) {
			locind = 4*ii + 2*jj + kk;
			vertind = block*PPB + 16*(bi+ii)+4*(bj+jj)+(bk+kk);
			mesh->connectivity[8*p+locind] = vertind;
		}

	}
#else
	// now build the mesh connectivity
	// trilinear elements naturally
	nside = floor(pow(mesh->npart+0.5, ONE_THIRD));
	for(i = 0; i < nside; ++i)
	for(j = 0; j < nside; ++j)
	for(k = 0; k < nside; ++k) {
		elemind = nside*nside*i + nside*j + k;
		for(ii = 0; ii < 2; ++ii)
		for(jj = 0; jj < 2; ++jj)
		for(kk = 0; kk < 2; ++kk) {
			locind = 4*ii + 2*jj + kk;
			vertind = nside*nside*((i+ii)%nside)
				+ nside*((j+jj)%nside) + ((k+kk)%nside);
			mesh->connectivity[8*elemind+locind] = vertind;
		}
	}
#endif

	psi_printf("...done.\n");
	return 1;
}




