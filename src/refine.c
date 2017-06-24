#include "refine.h"

// TODO: write an interpolation routing generally, for any order??
// Lagrange polynomials can be computed in log(n) time, recursively.
// Can the whole thing be done that way somehow?


// Indexing macros and connectivity for decomposing the unit cube into tetrahedra
#if PSI_NDIM == 3
#define ind2(ii, jj, kk) (4*(ii)+2*(jj)+(kk))
#define ind3(ii, jj, kk) (9*(ii)+3*(jj)+(kk))
#define ind5(ii, jj, kk) (25*(ii)+5*(jj)+(kk))
const static psi_int PSI_TET_INDS[6][4] = {{2, 0, 6, 18}, {8, 2, 6, 18}, 
		{8, 20, 2, 18}, {8, 24, 20, 18}, {8, 6, 24, 18}, {8, 26, 20, 24}};	
const static psi_int PSI_NTETS = 6;
#elif PSI_NDIM == 2
#define ind2(ii, jj) (2*(ii)+(jj))
#define ind3(ii, jj) (3*(ii)+(jj))
#define ind5(ii, jj) (5*(ii)+(jj))
const static psi_int PSI_TET_INDS[2][3] = {{0,2,6}, {2,6,8}};
const static psi_int PSI_NTETS = 2;
#endif
#define M_PER_TET (1.0/PSI_NTETS)

// private declarations
void psi_refine_triquadratic(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf, psi_int lvl);
void psi_refine_trilinear(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf, psi_int lvl);
void psi_interp_trilinear(psi_rvec* in, psi_rvec* out);
void psi_interp_triquadratic(psi_rvec* in, psi_rvec* out);
void psi_tets_from_3cube(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf);

void psi_tet_buffer_refine(psi_tet_buffer* tetbuf, psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_int order) {

	if(order == PSI_MESH_SIMPLEX) {
		// if it's a single tet, just copy it to the buffer as a dummy
		tetbuf->num = 1;
		memcpy(tetbuf->pos, pos, (PSI_NDIM+1)*sizeof(psi_rvec));
		memcpy(tetbuf->vel, vel, (PSI_NDIM+1)*sizeof(psi_rvec));
		tetbuf->mass[0] = mass;
	}
	else {
		// if it's a curvilinear element, refine to the buffer
		// TODO: make the interpolation code more homogeneous
		// TODO: when to use velocity??
		tetbuf->num = 0;
		if(order == PSI_MESH_LINEAR) {
			psi_rvec tmppos[27];
			psi_rvec tmpvel[27];
			psi_interp_trilinear(pos, tmppos);
			psi_interp_trilinear(vel, tmpvel);
			psi_refine_trilinear(tmppos, tmpvel, mass, tetbuf, 0);
		}
		else if(order == PSI_MESH_QUADRATIC) {
			psi_refine_triquadratic(pos, vel, mass, tetbuf, 0);
		}
	}
}

void psi_tet_buffer_init(psi_tet_buffer* buffer, psi_real tol, psi_int max_lvl) {
	// TODO: put a check on max_lvl??
	// TODO: heuristically make the intial capacity a bit smaller?
	buffer->num = 0;
	if(max_lvl > 6) {
		printf("Warning: for space, psi_tet_buffer.max_lvl is restricted to 6\n");
		max_lvl = 6;
	}
	buffer->max_lvl = max_lvl;
	buffer->ptol2 = tol*tol;
	buffer->capacity = (1 << PSI_NDIM*max_lvl)*PSI_NTETS; 
	buffer->pos = malloc(buffer->capacity*(PSI_NDIM+1)*sizeof(psi_rvec));
	buffer->vel = malloc(buffer->capacity*(PSI_NDIM+1)*sizeof(psi_rvec));
	buffer->mass = malloc(buffer->capacity*sizeof(psi_real));
	if(!buffer->pos || !buffer->vel || !buffer->mass) {
		printf("Failed to allocate psi_tet_buffer.\n");
	}
}

void psi_tet_buffer_destroy(psi_tet_buffer* buffer) {
	if(buffer->pos) free(buffer->pos);
	if(buffer->vel) free(buffer->vel);
	if(buffer->mass) free(buffer->mass);
	buffer->num = 0;
	buffer->capacity = 0;
}

void psi_tet_buffer_grow(psi_tet_buffer* buffer) {
	// growth factor of 1.5 seems to work well
	buffer->capacity = (psi_int) (1.5*buffer->capacity);
	psi_rvec* newpos = realloc(buffer->pos, buffer->capacity*(PSI_NDIM+1)*sizeof(psi_rvec));
	psi_rvec* newvel = realloc(buffer->vel, buffer->capacity*(PSI_NDIM+1)*sizeof(psi_rvec));
	psi_real* newmass = realloc(buffer->mass, buffer->capacity*sizeof(psi_real));
	if(newpos && newvel && newmass) {
		buffer->pos = newpos;
		buffer->vel = newvel;
		buffer->mass = newmass;
	}
	else {
		printf("Failed to grow psi_tet_buffer.\n");
	}
}

void psi_tets_from_3cube(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf) {

	psi_int tid, v, pind;
	psi_rvec* bufpos;
	psi_rvec* bufvel;

	// grow the buffer if needed
	if(tetbuf->num + PSI_NTETS >= tetbuf->capacity)
		psi_tet_buffer_grow(tetbuf);

	// push new tets onto the buffer
	for(tid = 0; tid < PSI_NTETS; ++tid) {
		bufpos = &tetbuf->pos[(PSI_NDIM+1)*tetbuf->num];
		bufvel = &tetbuf->vel[(PSI_NDIM+1)*tetbuf->num];
		for(v = 0; v < (PSI_NDIM+1); ++v) {
			pind = PSI_TET_INDS[tid][v];
			bufpos[v] = pos[pind];
			bufvel[v] = vel[pind];
		}
		tetbuf->mass[tetbuf->num] = mass*M_PER_TET;
		tetbuf->num++;
	}
}

void psi_refine_triquadratic(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf, psi_int lvl) {

	// TODO: these are different for the 48-tet connectivity!
#if PSI_NDIM == 3
	// all unique tet edges across the 3x3 unit cube. Order is first, midpoint, second.
	const static int edges3TQ[19][3] =  { { 0 , 1 , 2 }, { 0 , 3 , 6 }, { 0 , 9 , 18 }, { 2 , 4 , 6 }, 
											{ 2 , 5 , 8 }, { 2 , 10 , 18 }, { 2 , 11 , 20 }, { 6 , 7 , 8 },
										   	{ 6 , 12 , 18 }, { 6 , 15 , 24 }, { 8 , 13 , 18 }, { 8 , 14 , 20 }, 
											{ 8 , 16 , 24 }, { 8 , 17 , 26 }, { 18 , 19 , 20 }, { 18 , 21 , 24 }, 
											{ 20 , 22 , 24 }, { 20 , 23 , 26 }, { 24 , 25 , 26 } };
	const static int nedges3TQ = 19; 
#elif PSI_NDIM == 2
	const static psi_int edges3TQ[5][3] =  {{0,1,2}, {0,3,6}, {2,4,6}, {2,5,8}, {6,7,8}};
	const static psi_int nedges3TQ = 5; 
#endif

	psi_rvec interppos[125];
	psi_rvec interpvel[125];
	psi_rvec childpos[27];
	psi_rvec childvel[27];
	psi_int e, i, j, k, ii, jj, kk;
	psi_real midpt_tet, midpt_quad, err2;

	// check to see if all of the tet edges are within tolerance
	if(lvl < tetbuf->max_lvl) for(e = 0; e < nedges3TQ; ++e) {
		err2 = 0.0;
		for(i = 0; i < PSI_NDIM; ++i) {
			midpt_tet = 0.5*(pos[edges3TQ[e][0]].xyz[i] + pos[edges3TQ[e][2]].xyz[i]);
			midpt_quad = pos[edges3TQ[e][1]].xyz[i];
			err2 += (midpt_tet - midpt_quad)*(midpt_tet - midpt_quad);
		}

		// if this edge gives too large an error in position, stop and recurse
		if(err2 > tetbuf->ptol2) {

			// triquadratic interpolation
			psi_interp_triquadratic(pos, interppos);
			psi_interp_triquadratic(vel, interpvel); // TODO: option to use velocity

#if PSI_NDIM == 3
			// extract sub-cubes from interp
			for(i = 0; i < 2; ++i)
			for(j = 0; j < 2; ++j)
			for(k = 0; k < 2; ++k) {
				for(ii = 0; ii < 3; ++ii)
				for(jj = 0; jj < 3; ++jj)
				for(kk = 0; kk < 3; ++kk) {// TODO: option to use velocity
					childpos[ind3(ii,jj,kk)] = interppos[ind5(ii+2*i,jj+2*j,kk+2*k)];
					childvel[ind3(ii,jj,kk)] = interpvel[ind5(ii+2*i,jj+2*j,kk+2*k)];
				}
				psi_refine_triquadratic(childpos, childvel, 0.125*mass, tetbuf, lvl+1);
			}
#elif PSI_NDIM == 2
			for(i = 0; i < 2; ++i)
			for(j = 0; j < 2; ++j) {
				for(ii = 0; ii < 3; ++ii)
				for(jj = 0; jj < 3; ++jj) {
					childpos[ind3(ii,jj)] = interppos[ind5(ii+2*i,jj+2*j)];
					childvel[ind3(ii,jj)] = interpvel[ind5(ii+2*i,jj+2*j)];
				}
				psi_refine_triquadratic(childpos, childvel, 0.25*mass, tetbuf, lvl+1);
			}
#endif
			return;
		}
	}

	// if we have passed the tolerance test, refine the ucube to the buffer
	psi_tets_from_3cube(pos, vel, mass, tetbuf);
}

void psi_refine_trilinear(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_tet_buffer* tetbuf, psi_int lvl) {

	// TODO: these are different for the 48-tet connectivity!
#if PSI_NDIM == 3
	// all unique tet edges across the 3-cube. Order is first, midpoint, second.
	const static psi_int edges3TL[7][3] =  { { 2 , 4 , 6 }, { 2 , 10 , 18 }, 
										   	{ 6 , 12 , 18 }, { 8 , 13 , 18 }, { 8 , 14 , 20 }, 
											{ 8 , 16 , 24 }, { 20 , 22 , 24 } };
	const static psi_int nedges3TL = 7; 
#elif PSI_NDIM == 2
	const static psi_int edges3TL[5][3] =  {{0,1,2}, {0,3,6}, {2,4,6}, {2,5,8}, {6,7,8}};
	const static psi_int nedges3TL = 5; 
#endif

	psi_int e, i, j, k, ii, jj, kk;
	psi_rvec tmppos[27];
	psi_rvec tmpvel[27];
	psi_rvec childpos[8];
	psi_rvec childvel[8];
	psi_real midpt_tet, midpt_lin, err2;

	// check to see if all of the tet edges are within tolerance
	if(lvl < tetbuf->max_lvl) for(e = 0; e < nedges3TL; ++e) {
		err2 = 0.0;
		for(i = 0; i < PSI_NDIM; ++i) {
			midpt_tet = 0.5*(pos[edges3TL[e][0]].xyz[i] + pos[edges3TL[e][2]].xyz[i]);
			midpt_lin = pos[edges3TL[e][1]].xyz[i];
			err2 += (midpt_tet - midpt_lin)*(midpt_tet - midpt_lin);
		}

		// if this edge gives too large an error in position, stop and recurse
		if(err2 > tetbuf->ptol2) {

#if PSI_NDIM == 3
			for(i = 0; i < 2; ++i)
			for(j = 0; j < 2; ++j)
			for(k = 0; k < 2; ++k) {
				for(ii = 0; ii < 2; ++ii)
				for(jj = 0; jj < 2; ++jj)
				for(kk = 0; kk < 2; ++kk) {
					childpos[ind2(ii,jj,kk)] = pos[ind3(ii+i,jj+j,kk+k)];
					childvel[ind2(ii,jj,kk)] = vel[ind3(ii+i,jj+j,kk+k)];
				} 
				psi_interp_trilinear(childpos, tmppos);
				psi_interp_trilinear(childvel, tmpvel);
				psi_refine_trilinear(tmppos, tmpvel, 0.125*mass, tetbuf, lvl+1);
			}
#elif PSI_NDIM == 2
			for(i = 0; i < 2; ++i)
			for(j = 0; j < 2; ++j) {
				for(ii = 0; ii < 2; ++ii)
				for(jj = 0; jj < 2; ++jj) {
					childpos[ind2(ii,jj)] = pos[ind3(ii+i,jj+j)];
					childvel[ind2(ii,jj)] = vel[ind3(ii+i,jj+j)];
				}
				psi_interp_trilinear(childpos, tmppos);
				psi_interp_trilinear(childvel, tmpvel);
				psi_refine_trilinear(tmppos, tmpvel, 0.25*mass, tetbuf, lvl+1);
			}
#endif
			return;
		}
	}

	// if we have passed the tolerance test, refine to the buffer 
	psi_tets_from_3cube(pos, vel, mass, tetbuf);
}


void psi_interp_trilinear(psi_rvec* in, psi_rvec* out) {
	psi_int i, j, k, ax;
	for(ax = 0; ax < PSI_NDIM; ++ax) {
#if PSI_NDIM == 3
		for(i = 0; i < 2; ++i)	
		for(j = 0; j < 2; ++j)	
		for(k = 0; k < 2; ++k)	
			out[ind3(2*i, 2*j, 2*k)].xyz[ax] = in[ind2(i, j, k)].xyz[ax];
		for(j = 0; j < 3; j += 2)	
		for(k = 0; k < 3; k += 2)	
			out[ind3(1, j, k)].xyz[ax] = 0.5*(out[ind3(0, j, k)].xyz[ax] + out[ind3(2, j, k)].xyz[ax]);
		for(i = 0; i < 3; i += 1)	
		for(k = 0; k < 3; k += 2)	
			out[ind3(i, 1, k)].xyz[ax] = 0.5*(out[ind3(i, 0, k)].xyz[ax] + out[ind3(i, 2, k)].xyz[ax]);
		for(i = 0; i < 3; i += 1)	
		for(j = 0; j < 3; j += 1)	
			out[ind3(i, j, 1)].xyz[ax] = 0.5*(out[ind3(i, j, 0)].xyz[ax] + out[ind3(i, j, 2)].xyz[ax]);
#elif PSI_NDIM == 2
		for(i = 0; i < 2; ++i)	
		for(j = 0; j < 2; ++j)	
			out[ind3(2*i, 2*j)].xyz[ax] = in[ind2(i, j)].xyz[ax];
		for(j = 0; j < 3; j += 2)	
			out[ind3(1, j)].xyz[ax] = 0.5*(out[ind3(0, j)].xyz[ax] + out[ind3(2, j)].xyz[ax]);
		for(i = 0; i < 3; i += 1)	
			out[ind3(i, 1)].xyz[ax] = 0.5*(out[ind3(i, 0)].xyz[ax] + out[ind3(i, 2)].xyz[ax]);
#endif
	}
}

void psi_interp_triquadratic(psi_rvec* in, psi_rvec* out) {
	psi_int i, j, k, ax;
	for(ax = 0; ax < PSI_NDIM; ++ax) {
#if PSI_NDIM == 3
		for(i = 0; i < 3; ++i)	
		for(j = 0; j < 3; ++j)	
		for(k = 0; k < 3; ++k)	
			out[ind5(2*i, 2*j, 2*k)].xyz[ax] = in[ind3(i, j, k)].xyz[ax];
		for(j = 0; j < 5; j += 2)	
		for(k = 0; k < 5; k += 2) {
			out[ind5(1, j, k)].xyz[ax] = 0.125*(3.0*out[ind5(0, j, k)].xyz[ax] 
					+ 6.0*out[ind5(2, j, k)].xyz[ax] - out[ind5(4, j, k)].xyz[ax]);
			out[ind5(3, j, k)].xyz[ax] = 0.125*(3.0*out[ind5(4, j, k)].xyz[ax] 
					+ 6.0*out[ind5(2, j, k)].xyz[ax] - out[ind5(0, j, k)].xyz[ax]);
		}
		for(i = 0; i < 5; i += 1)	
		for(k = 0; k < 5; k += 2) {
			out[ind5(i, 1, k)].xyz[ax] = 0.125*(3.0*out[ind5(i, 0, k)].xyz[ax] 
					+ 6.0*out[ind5(i, 2, k)].xyz[ax] - out[ind5(i, 4, k)].xyz[ax]);
			out[ind5(i, 3, k)].xyz[ax] = 0.125*(3.0*out[ind5(i, 4, k)].xyz[ax] 
					+ 6.0*out[ind5(i, 2, k)].xyz[ax] - out[ind5(i, 0, k)].xyz[ax]);
		}
		for(i = 0; i < 5; i += 1)	
		for(j = 0; j < 5; j += 1) {
			out[ind5(i, j, 1)].xyz[ax] = 0.125*(3.0*out[ind5(i, j, 0)].xyz[ax] 
					+ 6.0*out[ind5(i, j, 2)].xyz[ax] - out[ind5(i, j, 4)].xyz[ax]);
			out[ind5(i, j, 3)].xyz[ax] = 0.125*(3.0*out[ind5(i, j, 4)].xyz[ax] 
					+ 6.0*out[ind5(i, j, 2)].xyz[ax] - out[ind5(i, j, 0)].xyz[ax]);
		}
#elif PSI_NDIM == 2
		for(i = 0; i < 3; ++i)	
		for(j = 0; j < 3; ++j)	
			out[ind5(2*i, 2*j)].xyz[ax] = in[ind3(i, j)].xyz[ax];
		for(j = 0; j < 5; j += 2) {	
			out[ind5(1, j)].xyz[ax] = 0.125*(3.0*out[ind5(0, j)].xyz[ax] 
					+ 6.0*out[ind5(2, j)].xyz[ax] - out[ind5(4, j)].xyz[ax]);
			out[ind5(3, j)].xyz[ax] = 0.125*(3.0*out[ind5(4, j)].xyz[ax] 
					+ 6.0*out[ind5(2, j)].xyz[ax] - out[ind5(0, j)].xyz[ax]);
		}
		for(i = 0; i < 5; i += 1) {	
			out[ind5(i, 1)].xyz[ax] = 0.125*(3.0*out[ind5(i, 0)].xyz[ax] 
					+ 6.0*out[ind5(i, 2)].xyz[ax] - out[ind5(i, 4)].xyz[ax]);
			out[ind5(i, 3)].xyz[ax] = 0.125*(3.0*out[ind5(i, 4)].xyz[ax] 
					+ 6.0*out[ind5(i, 2)].xyz[ax] - out[ind5(i, 0)].xyz[ax]);
		}
#endif
	}
}

#if 0 // 48-tet decomposition for testing
	const static psi_int conn[48][4] = {{13, 0, 9, 12}, {13, 0, 9, 10}, {13, 0, 3, 12}, {13, 0, 3, 4}, 
		{13, 0, 1, 10}, {13, 0, 1, 4}, {13, 2, 11, 14}, {13, 2, 11, 10}, {13, 2, 5, 14}, {13, 2, 5, 4}, 
		{13, 2, 1, 10}, {13, 2, 1, 4}, {13, 6, 15, 12}, {13, 6, 15, 16}, {13, 6, 3, 12}, {13, 6, 3, 4}, 
		{13, 6, 7, 16}, {13, 6, 7, 4}, {13, 8, 17, 14}, {13, 8, 17, 16}, {13, 8, 5, 14}, {13, 8, 5, 4}, 
		{13, 8, 7, 16}, {13, 8, 7, 4}, {13, 18, 9, 12}, {13, 18, 9, 10}, {13, 18, 21, 12}, {13, 18, 21, 22}, 
		{13, 18, 19, 10}, {13, 18, 19, 22}, {13, 20, 11, 14}, {13, 20, 11, 10}, {13, 20, 23, 14}, {13, 20, 23, 22}, 
		{13, 20, 19, 10}, {13, 20, 19, 22}, {13, 24, 15, 12}, {13, 24, 15, 16}, {13, 24, 21, 12}, {13, 24, 21, 22}, 
		{13, 24, 25, 16}, {13, 24, 25, 22}, {13, 26, 17, 14}, {13, 26, 17, 16}, {13, 26, 23, 14}, {13, 26, 23, 22}, 
		{13, 26, 25, 16}, {13, 26, 25, 22}};
	const static psi_int ntet = 48;
#endif

