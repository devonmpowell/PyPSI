/*
 *
 *		psi.c
 *		
 *		See psi.h for usage and documentation.
 *		
 *		Devon Powell
 *		11 Feb 2016
 *
 */		

#include "psi.h"

#include "rtree.h"
#include "grid.h"
#include "mesh.h"
#include "refine.h"
#include "geometry.h"
//#include "particles.h"

// internal top-level utilities
void psi_aabb(psi_rvec* pos, psi_int nverts, psi_rvec* rbox);
psi_int psi_aabb_periodic(psi_rvec* pos, psi_rvec* rbox, psi_rvec* window, psi_mesh* mesh); 
int psi_aabb_ixn(psi_rvec* rbox0, psi_rvec* rbox1, psi_rvec* ixn);
void psi_make_ghosts(psi_rvec* elems, psi_rvec* rboxes, psi_int* num, psi_int stride, psi_rvec* window, psi_mesh* mesh);

// psi implementation 
void psi_voxels(psi_grid* grid, psi_mesh* mesh, psi_rtree* rtree, psi_int mode, psi_real reftol, psi_int max_ref_lvl) {

	//setbuf(stdout, NULL);
	psi_int e, g, t, v, nghosts, tind, internal_rtree, e1, t1;

	// a local copy of the position in case it is modified due to periodicity
	psi_int vpere = mesh->elemtype;
	psi_rvec tpos[vpere], tvel[vpere], tpos1[vpere], tvel1[vpere];
	psi_real tmass, tmass1;
	psi_rvec trbox[2];
	psi_rvec *pos0, *vel0, *pos1, *vel1;
	psi_rvec rbox0[2], rbox1[2], mbox[2];
	psi_real mass0, mass1;

	// ghost tetrahedra for handling periodicity
	psi_rvec gpos[(1<<mesh->dim)*(mesh->dim+1)];
	psi_rvec grbox[(1<<mesh->dim)*2];

	// constant-size buffer for max refinement level
	psi_rtree myrtree;
	psi_rtree_query qry;
	psi_tet_buffer tetbuf, tetbuf1;
	psi_tet_buffer_init(&tetbuf, reftol, max_ref_lvl);
	if(mode == PSI_MODE_ANNIHILATION)
		psi_tet_buffer_init(&tetbuf1, reftol, max_ref_lvl);

	// check to see if we must build an rtree for this function
	// if so, create one from the mesh
	internal_rtree = (mode == PSI_MODE_ANNIHILATION && rtree == NULL);
	if(internal_rtree) {
		psi_printf("No R*-tree given. Building one now.\n");
		rtree = &myrtree; 
		psi_rtree_from_mesh(rtree, mesh, grid->window);
	}

	// loop over all elements
	for(e = 0; e < mesh->nelem; ++e) {

		// a local copy of the element, in case it's modified
		// make it periodic, get its bounding box, and check it against the grid
		tmass = mesh->mass[e]; 
		for(v = 0; v < vpere; ++v) {
			tind = mesh->connectivity[e*vpere+v];
			tpos[v] = mesh->pos[tind];
			tvel[v] = mesh->vel[tind];
		}
		if(!psi_aabb_periodic(tpos, trbox, grid->window, mesh)) 
			continue;

		// refine elements into the tet buffer 
		psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, mesh->elemtype);

		// linear in the density
		if(mode == PSI_MODE_DENSITY) {

			// loop over each tet in the buffer
			// copy the position to the ghost array and compute its aabb
			// make ghosts and sample tets to the grid
			for(t = 0; t < tetbuf.num; ++t) {
				memcpy(gpos, &tetbuf.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
				psi_aabb(gpos, mesh->dim+1,grbox);
				psi_make_ghosts(gpos, grbox, &nghosts, (mesh->dim+1), grid->window, mesh);
				for(g = 0; g < nghosts; ++g) 
					psi_voxelize_tet(&gpos[(mesh->dim+1)*g], &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
			}
		}

		// quadratic in the density
		else if(mode == PSI_MODE_ANNIHILATION) {

			// TODO: make and save the clip planes, density, rboxes for tetbuf 1 here

			// query the rtree for overlapping elements,
			// refine them, and intersect all pairs of resulting tets
			// TODO: ghosts and periodicity
			psi_rtree_query_init(&qry, rtree, trbox);
			while(psi_rtree_query_next(&qry, &e1)) {

				// get the overlapping element and refine it
				tmass1 = mesh->mass[e1]; 
				for(v = 0; v < vpere; ++v) {
					tind = mesh->connectivity[e1*vpere+v];
					tpos1[v] = mesh->pos[tind];
					tvel1[v] = mesh->vel[tind];
				}
				psi_tet_buffer_refine(&tetbuf1, tpos1, tvel1, tmass1, mesh->elemtype);

				// loop over each pair of tets between the refine buffers 
				for(t = 0; t < tetbuf.num; ++t) {

					// get the bounding box of the first tet
					// to use for the rtree query to find overlapping tets
					pos0 = &tetbuf.pos[(PSI_NDIM+1)*t];
					vel0 = &tetbuf.vel[(PSI_NDIM+1)*t];
					mass0 = tetbuf.mass[t];
					psi_aabb(pos0, PSI_NDIM+1, rbox0); 

					// loop over the second buffer
					for(t1 = 0; t1 < tetbuf1.num; ++t1) {
						pos1 = &tetbuf1.pos[(PSI_NDIM+1)*t1];
						vel1 = &tetbuf1.vel[(PSI_NDIM+1)*t1];
						mass1 = tetbuf1.mass[t1];
						psi_aabb(pos1, PSI_NDIM+1, rbox1); 

						// voxelize the annihilation rate between pairs of tets
						// only of mbox (the mutual intersection of bounding boxes) exists
						if(psi_aabb_ixn(rbox0, rbox1, mbox))
							psi_voxelize_annihilation(pos0, vel0, mass0, pos1, vel1, mass1, mbox, grid);
					}
				}
			}
		}

		if(e%PRINT_EVERY==0)
			psi_printf("\rElement %d of %d, %.1f%%", e, mesh->nelem, (100.0*e)/mesh->nelem);
#ifdef PYMODULE
		PyErr_CheckSignals(); // TODO: why does this do nothing??
#endif
	}
	psi_printf("\n");

	// free temporary buffers
	psi_tet_buffer_destroy(&tetbuf);
	if(mode == PSI_MODE_ANNIHILATION)
		psi_tet_buffer_destroy(&tetbuf1);
	if(internal_rtree) 
		psi_rtree_destroy(rtree);

}

void psi_sample_vdf(psi_rvec samppos, psi_mesh* mesh, psi_rtree* rtree, 
		psi_real** rhoout, psi_rvec** velout, psi_int* nsamp , psi_real reftol, psi_int max_ref_lvl) {

	psi_int i, j, t, e, v, tind, internal_rtree, capacity;
	psi_rvec *pos, *vel;
	psi_real mass, vol;
	psi_real bcoords[(PSI_NDIM+1)];

	psi_int vpere = mesh->elemtype;
	psi_rvec tpos[vpere], tvel[vpere];
	psi_real tmass;

	// set up the return buffers
	// NOTE: It is up to the user or to Python
	// to manage these buffers after they are returned!!!
	capacity = 1024; 
	*rhoout = (psi_real*) psi_malloc(capacity*sizeof(psi_real));
	*velout = (psi_rvec*) psi_malloc(capacity*sizeof(psi_rvec));
	*nsamp = 0;

	// constant-size buffer for max refinement level
	psi_tet_buffer tetbuf;
	psi_tet_buffer_init(&tetbuf, reftol, max_ref_lvl);

	// check to see if we must build an rtree for this function
	// if so, create one from the mesh
	psi_rtree myrtree;
	internal_rtree = (rtree == NULL);
	if(internal_rtree) {
		psi_printf("No R*-tree given. Building one now.\n");
		rtree = &myrtree; 
		psi_rtree_from_mesh(rtree, mesh, mesh->box);
	}

	// query the r*-tree using a degenerate query box
	psi_rtree_query qry;
	psi_rvec qbox[2];
	qbox[0] = samppos;
	qbox[1] = samppos;
	psi_rtree_query_init(&qry, rtree, qbox);
	while(psi_rtree_query_next(&qry, &e)) {

		// get the overlapping element 
		// refine the element into tetrahedra
		tmass = mesh->mass[e]; 
		for(v = 0; v < vpere; ++v) {
			tind = mesh->connectivity[e*vpere+v];
			tpos[v] = mesh->pos[tind];
			tvel[v] = mesh->vel[tind];
		}
		psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, mesh->elemtype);

		// loop over each tet in the buffer and sample the vdf
		for(t = 0; t < tetbuf.num; ++t) {
			pos = &tetbuf.pos[(PSI_NDIM+1)*t];
			vel = &tetbuf.vel[(PSI_NDIM+1)*t];
			mass = tetbuf.mass[t];

			// get the barycentric coordinates from the matrix inverse
			// return if the point lies outside any face of the tet
			// interpolate velocity and add the sample to the buffer
			vol = psi_barycentric(pos, samppos, bcoords);

			// skip fully degenerate tets
			// skip if the sample point lies outside of the tet
			if(vol <= 0.0) goto next_tet; 
			for(i = 0; i < PSI_NDIM+1; ++i) 
				if(bcoords[i] < 0.0) goto next_tet; 

			// grow the buffer if needed
			if(*nsamp >= capacity) {
				capacity *= 2; 
				*rhoout = (psi_real*) psi_realloc(*rhoout, capacity*sizeof(psi_real));
				*velout = (psi_rvec*) psi_realloc(*velout, capacity*sizeof(psi_rvec));
			}

			// use the barycentric coordinates to sample the vdf and get the density
			for(i = 0; i < PSI_NDIM; ++i) {
				(*velout)[*nsamp].xyz[i] = 0.0; 
				for(j = 0; j < PSI_NDIM+1; ++j)
					(*velout)[*nsamp].xyz[i] += bcoords[j]*vel[j].xyz[i];
			}
			(*rhoout)[*nsamp] = mass/vol; // vol is guaranteed positive from psi_barycentric
			(*nsamp)++;
			next_tet: continue;
		}
	}

	// free stuff
	psi_tet_buffer_destroy(&tetbuf);
	if(internal_rtree) 
		psi_rtree_destroy(rtree);
}

void psi_aabb(psi_rvec* pos, psi_int nverts, psi_rvec* rbox) {
	// get the bounding box 
	psi_int i, v;
	for(i = 0; i < 3; ++i) {
		rbox[0].xyz[i] = 1.0e30;
		rbox[1].xyz[i] = -1.0e30;
	}
	for(v = 0; v < nverts; ++v)
	for(i = 0; i < 3; ++i) {
		// TODO: This is naive wrt quadratic elements!
		// Could cause mass loss around the edges
		if(pos[v].xyz[i] < rbox[0].xyz[i]) rbox[0].xyz[i] = pos[v].xyz[i];
		if(pos[v].xyz[i] > rbox[1].xyz[i]) rbox[1].xyz[i] = pos[v].xyz[i];
	}
}

int psi_aabb_ixn(psi_rvec* rbox0, psi_rvec* rbox1, psi_rvec* ixn) {

	// constructs the intersection bounding box
	// returns 1 if they intersect, 0 if not
	int i, ret;
	for(i = 0; i < 3; ++i) {
		ixn[0].xyz[i] = fmax(rbox0[0].xyz[i], rbox1[0].xyz[i]);
		ixn[1].xyz[i] = fmin(rbox0[1].xyz[i], rbox1[1].xyz[i]);
	}
	for(i = 0, ret = 1; i < 3; ++i) {
		ret &= (ixn[1].xyz[i] > ixn[0].xyz[i]);
	}
	return ret;
}


psi_int psi_aabb_periodic(psi_rvec* pos, psi_rvec* rbox, psi_rvec* window, psi_mesh* mesh) {

	// pos and rbox can both be modified!
	psi_int i, v;
	psi_real span;

	psi_aabb(pos, mesh->elemtype, rbox);

	// if the grid is periodic, move the vertices as needed
	if(mesh->periodic) {
		for(i = 0; i < mesh->dim; ++i) {
			span = mesh->box[1].xyz[i]-mesh->box[0].xyz[i];
			if(rbox[1].xyz[i]-rbox[0].xyz[i] > 0.5*span) {
				rbox[0].xyz[i] = 1.0e30;
				rbox[1].xyz[i] = -1.0e30;
				for(v = 0; v < mesh->elemtype; ++v) {
					// wrap elements that cross the right-hand boundary while recomputing the box
					if(pos[v].xyz[i] > mesh->box[0].xyz[i] + 0.5*span) pos[v].xyz[i] -= span;
					if(pos[v].xyz[i] < rbox[0].xyz[i]) rbox[0].xyz[i] = pos[v].xyz[i];
					if(pos[v].xyz[i] > rbox[1].xyz[i]) rbox[1].xyz[i] = pos[v].xyz[i];
				}
			}
			if(rbox[0].xyz[i] > window[1].xyz[i]) return 0; 
			if(rbox[1].xyz[i] < window[0].xyz[i]
				&& rbox[0].xyz[i]+span > window[1].xyz[i]) return 0; 
		}
	}
	else for(i = 0; i < mesh->dim; ++i) {
		if(rbox[0].xyz[i] > window[1].xyz[i]) return 0; 
		if(rbox[1].xyz[i] < window[0].xyz[i]) return 0;
	}
	return 1;
}

void psi_make_ghosts(psi_rvec* elems, psi_rvec* rboxes, psi_int* num, psi_int stride, psi_rvec* window, psi_mesh* mesh) {

	psi_int i, p, v, pflags;
	psi_real pshift;
	psi_rvec span;

	// if periodic, make ghosts on the fly
	(*num) = 1;
	pflags = 0;
	if(mesh->periodic) {
		for(i = 0; i < mesh->dim; ++i) {
			span.xyz[i] = mesh->box[1].xyz[i]-mesh->box[0].xyz[i]; 
			if(rboxes[0].xyz[i] < mesh->box[0].xyz[i]) pflags |= (1 << i);
		}
		for(p = 1; p < (1<<mesh->dim); ++p) {
			if((pflags & p) == p) {
				for(i = 0; i < mesh->dim; ++i) {
					pshift = span.xyz[i]*(1 & (p >> i));
					rboxes[2*(*num)+0].xyz[i] = rboxes[0].xyz[i] + pshift; 
					rboxes[2*(*num)+1].xyz[i] = rboxes[1].xyz[i] + pshift;
					if(rboxes[2*(*num)+0].xyz[i] > window[1].xyz[i]) goto next_shift;
					if(rboxes[2*(*num)+1].xyz[i] < window[0].xyz[i]) goto next_shift;
					for(v = 0; v < stride; ++v)
						elems[stride*(*num)+v].xyz[i] = elems[v].xyz[i] + pshift; 
				}
				(*num)++;
			}
			next_shift: continue;
		}	
	}
}


#if 0

void psi_count_streams(psi_rtree* rtree, psi_rvec* samppos, psi_int* count, psi_real* rho, psi_int nsamp, psi_real reftol) {

	// NOTE: mass, vel are dummies here
	psi_int i, j, t, s, e;
	psi_rvec *pos, vel[27];
	psi_real mass, vol;
	psi_real bcoords[(PSI_NDIM+1)];

	// constant-size buffer for max refinement level
	// TODO: user-issued max level
	psi_tet_buffer tetbuf;
	psi_int max_lvl = 4;
	if(rtree->order == 0) max_lvl = 0;
	psi_tet_buffer_init(&tetbuf, reftol, max_lvl);

	psi_rtree_query qry;
	psi_rvec qbox[2];
	for(s = 0; s < nsamp; ++s) {
		psi_printf("Sample %d of %d (%.1f\%)\r", s, nsamp, 100.0*s/nsamp);
		
		// query the r*-tree using a degenerate query box
		qbox[0] = samppos[s];
		qbox[1] = samppos[s];
		psi_rtree_query_init(&qry, rtree, qbox);
		while(psi_rtree_query_next(&qry, &pos, &vel, &mass, &e)) {
	
#if 1
			// refine the element into tetrahedra
			psi_tet_buffer_refine(&tetbuf, pos, vel, mass, rtree->order);
	
			// loop over each tet in the buffer and sample the vdf
			for(t = 0; t < tetbuf.num; ++t) {
				pos = &tetbuf.pos[(PSI_NDIM+1)*t];
				mass = tetbuf.mass[t];
	
				// get the barycentric coordinates from the matrix inverse
				// return if the point lies outside any face of the tet
				// interpolate velocity and add the sample to the buffer
				vol = psi_barycentric(pos, samppos[s], bcoords);
	
				// skip fully degenerate tets
				// skip if the sample point lies outside of the tet
				if(!vol) goto next_tet; 
				for(i = 0; i < PSI_NDIM+1; ++i) if(bcoords[i] < 0.0) goto next_tet;

				// increment the counter if we intersect a stream
				count[s]++;
				rho[s] += mass/vol; // vol is guaranteed positive from psi_barycentric
				next_tet: continue;
			}
#else
			count[s]++;
#endif
		}
	}
	psi_printf("\n");
	psi_tet_buffer_destroy(&tetbuf);
}


void psi_particle_mesh(psi_rvec* pos, psi_rvec* vel, psi_real* mass, psi_int npart, psi_dest_grid* grid) {

	// TODO: add other particle cloud types?

	psi_int t, i, p, pflags;
	psi_real pshift;
	psi_rvec span;
	psi_rvec tpos, ppos;

	// compute the box span
	if(grid->periodic) for(i = 0; i < PSI_NDIM; ++i)
		span.xyz[i] = grid->box[1][i]-grid->box[0][i]; 

	// loop over all elements
	for(t = 0; t < npart; ++t) {

		// get a local copy 
		tpos = pos[t];

		// check for periodicity and make ghosts if needed 
		if(grid->periodic) {
			pflags = 0;
			for(i = 0; i < PSI_NDIM; ++i) {
				if(tpos.xyz[i] + 0.5*grid->d[i] >= grid->box[1][i])
					tpos.xyz[i] -= span.xyz[i];
				if(tpos.xyz[i] - 0.5*grid->d[i] < grid->box[0][i]) 
					pflags |= (1 << i);
			}
			if(pflags) for(p = 1; p < (1<<PSI_NDIM); ++p) {
				if((pflags & p) == p) {
					for(i = 0; i < PSI_NDIM; ++i) {
						pshift = span.xyz[i]*(1 & (p >> i));
						ppos.xyz[i] = tpos.xyz[i] + pshift;
					}
					psi_cic(ppos, vel[t], mass[t], grid);
				}
			}	
		}
		psi_cic(tpos, vel[t], mass[t], grid);
	}
}

void psi_element_block_from_grid(psi_rvec* pos_in, psi_rvec* vel_in, psi_real mpere, psi_dvec* ngridp, 
		psi_int order, psi_dvec* startp, psi_dvec* nump, psi_rvec* pos_out, psi_rvec* vel_out, psi_real* mass_out) {

	psi_int i, j, k, ii, jj, kk, e, v, vpere;
	psi_rvec* curpos;
	psi_rvec* curvel;

	// get the actual vectors
	// passing from Python needs pointers
	psi_dvec ngrid = *ngridp;
	psi_dvec start = *startp;
	psi_dvec num = *nump;

	vpere = PSI_VERTS_PER_ELEM[order];


	// TODO: calculate strides first!
	
#if PSI_NDIM == 3
#define inind(ii, jj, kk) (ngrid.j*ngrid.k*(ii)+ngrid.k*(jj)+(kk))
	for(i = 0, e = 0; i < num.i; ++i)
	for(j = 0; j < num.j; ++j)
	for(k = 0; k < num.k; ++k, ++e) {
		curpos = &pos_out[e*vpere];
		curvel = &vel_out[e*vpere];
		mass_out[e] = mpere;
		for(ii = 0; ii < order+1; ++ii)
		for(jj = 0; jj < order+1; ++jj)
		for(kk = 0; kk < order+1; ++kk) {
			v = (order+1)*(order+1)*ii+(order+1)*jj+kk;
			curpos[v] = pos_in[inind(((i+start.i)*order+ii)%ngrid.i, 
					((j+start.j)*order+jj)%ngrid.j, ((k+start.k)*order+kk)%ngrid.k)];
			curvel[v] = vel_in[inind(((i+start.i)*order+ii)%ngrid.i, 
					((j+start.j)*order+jj)%ngrid.j, ((k+start.k)*order+kk)%ngrid.k)];
		}
	}	
#elif PSI_NDIM == 2
#define inind(ii, jj) (ngrid.j*(ii)+(jj))
	for(i = 0, e = 0; i < num.i; ++i)
	for(j = 0; j < num.j; ++j, ++e) {
		curpos = &pos_out[e*vpere];
		curvel = &vel_out[e*vpere];
		mass_out[e] = mpere;
		for(ii = 0; ii < order+1; ++ii)
		for(jj = 0; jj < order+1; ++jj) {
			v = (order+1)*ii+jj;
			curpos[v] = pos_in[inind(((i+start.i)*order+ii)%ngrid.i, 
					((j+start.j)*order+jj)%ngrid.j)];
			curvel[v] = vel_in[inind(((i+start.i)*order+ii)%ngrid.i, 
					((j+start.j)*order+jj)%ngrid.j)];
		}
	}	
#endif

}

#if PSI_NDIM == 3
void psi_element_blocks_from_lg3(float* pos_in, float* vel_in, float* mass_in, psi_int npatch, 
		psi_int order, psi_rvec* pos_out, psi_rvec* vel_out, psi_real* mass_out) {

	psi_int c, p, i, vpere, ost, epbk;
	vpere = PSI_VERTS_PER_ELEM[order];

	// for passing to element_blocks_from_grid()
	psi_rvec tmppos[125], tmpvel[125];
	psi_real tmpmass;
	psi_dvec elemstart, nelem, ngrid;
	for(c = 0; c < 3; ++c) {
		nelem.ijk[c] = 4/order;
		elemstart.ijk[c] = 0;
		ngrid.ijk[c] = 5;
	}
	epbk = nelem.i*nelem.j*nelem.k;

	// loop over 5x5x5 Lagrangian patches
	for(p = 0; p < npatch; ++p) {
		// copy to tmp arrays, converting to proper data type if necessary
		for(i = 0; i < 125; ++i)
		for(c = 0; c < 3; ++c) {
			tmppos[i].xyz[c] = pos_in[3*128*p + 3*i + c];
			tmpvel[i].xyz[c] = vel_in[3*128*p + 3*i + c];
		}
		tmpmass = mass_in[p];
		ost = p*epbk;
		psi_element_block_from_grid(tmppos, tmpvel, tmpmass/epbk,
				&ngrid, order, &elemstart, &nelem, &pos_out[ost*vpere], &vel_out[ost*vpere], &mass_out[ost]);
	}
}
#endif

#endif


