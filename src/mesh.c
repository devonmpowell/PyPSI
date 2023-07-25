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

	unsigned int blksize;
	int e, i;
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

	unsigned int blksize;
	int p, e, nside, ii, jj, kk, i, j, k;
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


static const psi_int gevolution_ntile = 4;
static const psi_int gevolution_ntile3 = 64;
static const psi_int gevolution_flat_to_cube[64][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}, {1, 0, 0}, {0, 0, 2}, {0, 0, 3}, {0, 1, 3}, {0, 1, 2}, {1, 1, 2}, {1, 1, 3}, {1, 0, 3}, {1, 0, 2}, {0, 2, 2}, {0, 2, 3}, {0, 3, 3}, {0, 3, 2}, {1, 3, 2}, {1, 3, 3}, {1, 2, 3}, {1, 2, 2}, {0, 2, 0}, {0, 2, 1}, {0, 3, 1}, {0, 3, 0}, {1, 3, 0}, {1, 3, 1}, {1, 2, 1}, {1, 2, 0}, {2, 2, 0}, {2, 2, 1}, {2, 3, 1}, {2, 3, 0}, {3, 3, 0}, {3, 3, 1}, {3, 2, 1}, {3, 2, 0}, {2, 2, 2}, {2, 2, 3}, {2, 3, 3}, {2, 3, 2}, {3, 3, 2}, {3, 3, 3}, {3, 2, 3}, {3, 2, 2}, {2, 0, 2}, {2, 0, 3}, {2, 1, 3}, {2, 1, 2}, {3, 1, 2}, {3, 1, 3}, {3, 0, 3}, {3, 0, 2}, {2, 0, 0}, {2, 0, 1}, {2, 1, 1}, {2, 1, 0}, {3, 1, 0}, {3, 1, 1}, {3, 0, 1}, {3, 0, 0}};
static const psi_int gevolution_cube_to_flat[4][4][4] = {{{ 0, 1, 8, 9 }, { 3, 2, 11, 10 }, { 24, 25, 16, 17 }, { 27, 26, 19, 18 } }, { { 7, 6, 15, 14 }, { 4, 5, 12, 13 }, { 31, 30, 23, 22 }, { 28, 29, 20, 21 } }, { { 56, 57, 48, 49 }, { 59, 58, 51, 50 }, { 32, 33, 40, 41 }, { 35, 34, 43, 42 } }, { { 63, 62, 55, 54 }, { 60, 61, 52, 53 }, { 39, 38, 47, 46 }, { 36, 37, 44, 45 } } };

void gevolution_lagrangian_cube_indices(int global_id, int tile_nside, int* cube_ids) {

	int tile_ns2, tile_idx, part_idx, flat_tile_idx;
	int ti, tj, tk, pi, pj, pk, pii, pjj, pkk, dti, dtj, dtk, ii, jj, kk; 

	// generate the mesh connectivity
	tile_ns2 = tile_nside*tile_nside; 

	// index of the tile
	tile_idx = global_id/gevolution_ntile3;
	tk = tile_idx/tile_ns2;
	tj = (tile_idx - tile_ns2*tk)/tile_nside;
	ti = (tile_idx - tile_ns2*tk - tile_nside*tj);

	// particle indices within the tile
	part_idx = global_id%gevolution_ntile3;
	pi = gevolution_flat_to_cube[part_idx][0]; 
	pj = gevolution_flat_to_cube[part_idx][1]; 
	pk = gevolution_flat_to_cube[part_idx][2]; 

	// build a logical cube of 8 neighbor particles
	for(ii = 0; ii < 2; ++ii)
	for(jj = 0; jj < 2; ++jj)
	for(kk = 0; kk < 2; ++kk) {

		// get the neighbor particle/block in Cartesian indices
		pii = pi+ii;
		pjj = pj+jj;
		pkk = pk+kk;
		dti = pii/gevolution_ntile; 
		dtj = pjj/gevolution_ntile;
		dtk = pkk/gevolution_ntile;
		pii %= gevolution_ntile;
		pjj %= gevolution_ntile;
		pkk %= gevolution_ntile;
		
		// the tile index and global particle index
		flat_tile_idx = ((ti+dti)%tile_nside) + tile_nside*((tj+dtj)%tile_nside) + tile_ns2*((tk+dtk)%tile_nside);
		cube_ids[4*ii + 2*jj + kk] = flat_tile_idx*gevolution_ntile3 + gevolution_cube_to_flat[pii][pjj][pkk];
	}
}




int load_gevolution(psi_mesh* mesh, const char* filename) {

	unsigned int blksize; 
	int p, e, i, nside;
	int tile_nside;
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
	printf("z = %f\n", header.redshift);
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

	// generate the mesh connectivity
	nside = floor(pow(mesh->npart+0.5, ONE_THIRD));
	tile_nside = nside/gevolution_ntile;
	for(p = 0; p < mesh->npart; ++p) {
		gevolution_lagrangian_cube_indices(p, tile_nside, &mesh->connectivity[8*p]);
	}

	psi_printf("...done.\n");
	return 1;
}




