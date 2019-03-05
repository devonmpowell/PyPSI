#include "geometry.h"

// geometry-related constants
#if PSI_NDIM == 3
#define PSI_TET_FAC ONE_SIXTH 
#elif PSI_NDIM == 2
#define PSI_TET_FAC 0.5
#endif

#define min(x, y) (((x) < (y))? (x) : (y))
#define max(x, y) (((x) > (y))? (x) : (y))

// struct to voxelize a polyhedron
// contains a stack and some grid information
typedef struct {
	psi_poly stack[64];
	psi_int nstack;
	psi_int polyorder;
	psi_grid* grid;
} psi_voxels;
void psi_voxels_init(psi_voxels* vox, psi_poly* poly, psi_int polyorder, psi_rvec* rbox, psi_grid* grid);
psi_int psi_voxels_next(psi_voxels* vox, psi_real* moments, psi_int* gridind);

// internal declarations for very low-level voxelization routines
void psi_clip(psi_poly* poly, psi_plane* planes, psi_int nplanes);
void psi_split_coord(psi_poly* inpoly, psi_poly* outpolys, psi_real coord, psi_int ax);
void psi_reduce(psi_poly* poly, psi_real* moments, psi_int polyorder, psi_int weight);
void psi_init_tet(psi_poly* poly, psi_rvec* verts);
psi_real psi_orient_tet(psi_rvec* pos, psi_rvec* vel);


// the top-level voxelization routine
void psi_voxelize_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_grid* grid) {

	psi_int i, j, ii, jj, gridind, polyorder;
	psi_real cm[PSI_NDIM], cv[PSI_NDIM], sxx[PSI_NDIM][PSI_NDIM], sxv[PSI_NDIM][PSI_NDIM], svv[PSI_NDIM][PSI_NDIM];
	psi_real cwght, dwght;
	psi_real DX[PSI_NDIM][PSI_NDIM], DV[PSI_NDIM][PSI_NDIM];
	psi_real moments[10];
	psi_poly curpoly;

#if PSI_NDIM == 3
	const static psi_int qinds[3][3] = {{4,5,6},{5,7,8},{6,8,9}};
#elif PSI_NDIM == 2
	const static psi_int qinds[2][2] = {{3,4},{4,5}};
#endif

	// set mass to the volume if volume-weighted
	// TODO: this swaps vertices if needed. Maybe do a more naive implementation?
	//if(grid->weight == PSI_GRID_MEIGHT_VOL) mass = psi_orient_tet(pos, vel); 

	// deformation matrices for the input tet's mass coordinates
	for(i = 0; i < PSI_NDIM; ++i)
	for(j = 0; j < PSI_NDIM; ++j) {
		DX[i][j] = pos[j+1].xyz[i]-pos[0].xyz[i];
		DV[i][j] = vel[j+1].xyz[i]-vel[0].xyz[i];
	}

	// initialize the tet as an edge-vertex graph 
	// and clamp it to the grid
	psi_init_tet(&curpoly, pos);

	// determine the correct poly order
	polyorder = 0;
	if(grid->fields[PSI_GRID_X] || grid->fields[PSI_GRID_V]) {
		polyorder = 1;	
	}
	if(grid->fields[PSI_GRID_XX] || grid->fields[PSI_GRID_XV] || grid->fields[PSI_GRID_VV]) {
		polyorder = 2;	
	}

	// initialize the voxelization iterator
	psi_voxels vox;
	psi_voxels_init(&vox, &curpoly, polyorder, rbox, grid);
	while(psi_voxels_next(&vox, moments, &gridind)) {
	
		// get moments for this fragment from the mass coordinates 
		if(grid->fields[PSI_GRID_X])
			for(i = 0; i < PSI_NDIM; ++i) {
				cm[i] = moments[0]*pos[0].xyz[i];
				for(j = 0; j < PSI_NDIM; ++j)
					cm[i] += moments[j+1]*DX[i][j];
				cm[i] /= moments[0];
			}
		if(grid->fields[PSI_GRID_V])
			for(i = 0; i < PSI_NDIM; ++i) {
				cv[i] = moments[0]*vel[0].xyz[i];
				for(j = 0; j < PSI_NDIM; ++j) 
					cv[i] += moments[j+1]*DV[i][j];
				cv[i] /= moments[0];
			}
		if(grid->fields[PSI_GRID_XX]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				sxx[i][j] = moments[0]*(pos[0].xyz[i]-cm[i])*(pos[0].xyz[j]-cm[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					sxx[i][j] += DX[i][ii]*moments[1+ii]*(pos[0].xyz[j]-cm[j]);
					sxx[i][j] += DX[j][ii]*moments[1+ii]*(pos[0].xyz[i]-cm[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					sxx[i][j] += DX[i][ii]*DX[j][jj]*moments[qinds[ii][jj]];
				sxx[i][j] /= moments[0];
			}
		}
		if(grid->fields[PSI_GRID_XV]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				sxv[i][j] = moments[0]*(pos[0].xyz[i]-cm[i])*(vel[0].xyz[j]-cv[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					sxv[i][j] += DX[i][ii]*moments[1+ii]*(pos[0].xyz[j]-cm[j]);
					sxv[i][j] += DV[j][ii]*moments[1+ii]*(vel[0].xyz[i]-cv[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					sxv[i][j] += DX[i][ii]*DV[j][jj]*moments[qinds[ii][jj]];
				sxv[i][j] /= moments[0];
			}
		}
		if(grid->fields[PSI_GRID_VV]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				svv[i][j] = moments[0]*(vel[0].xyz[i]-cv[i])*(vel[0].xyz[j]-cv[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					svv[i][j] += DV[i][ii]*moments[1+ii]*(vel[0].xyz[j]-cv[j]);
					svv[i][j] += DV[j][ii]*moments[1+ii]*(vel[0].xyz[i]-cv[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					svv[i][j] += DV[i][ii]*DV[j][jj]*moments[qinds[ii][jj]];
				svv[i][j] /= moments[0];
			}
		}

		// update everything online 
		// TODO: go back to adding up just the raw moments...
		cwght = mass*moments[0]/(grid->fields[PSI_GRID_M][gridind] + mass*moments[0]); 
		dwght = grid->fields[PSI_GRID_M][gridind]/(grid->fields[PSI_GRID_M][gridind] + mass*moments[0]); 
		if(grid->fields[PSI_GRID_XX]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_XX][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_XX][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*sxx[i][j] 
				+ cwght*dwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i])*(cm[j]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_XV]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_XV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_XV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*sxv[i][j] 
				+ cwght*dwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i])*(cv[j]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_VV]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_VV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_VV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*svv[i][j] 
				+ cwght*dwght*(cv[i]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i])*(cv[j]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_X]) for(i = 0; i < PSI_NDIM; ++i) 
			grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i] += cwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i]);
		if(grid->fields[PSI_GRID_V]) for(i = 0; i < PSI_NDIM; ++i) 
			grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i] += cwght*(cv[i]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i]);

		if(grid->fields[PSI_GRID_M]) grid->fields[PSI_GRID_M][gridind] += mass*moments[0];
	}
}

psi_real psi_orient_tet(psi_rvec* pos, psi_rvec* vel) {
	psi_real adx, bdx, cdx;
	psi_real ady, bdy, cdy;
	psi_real adz, bdz, cdz;
	psi_real vol;
	psi_rvec swp;
#if 1 
	adx = pos[0].x - pos[3].x;
	bdx = pos[1].x - pos[3].x;
	cdx = pos[2].x - pos[3].x;
	ady = pos[0].y - pos[3].y;
	bdy = pos[1].y - pos[3].y;
	cdy = pos[2].y - pos[3].y;
	adz = pos[0].z - pos[3].z;
	bdz = pos[1].z - pos[3].z;
	cdz = pos[2].z - pos[3].z;
	vol = -ONE_SIXTH*(adx*(bdy*cdz-bdz*cdy)
			+bdx*(cdy*adz-cdz*ady)+cdx*(ady*bdz-adz*bdy));
#elif 0 
	adx = pos[0].x - pos[2].x;
	bdx = pos[1].x - pos[2].x;
	ady = pos[0].y - pos[2].y;
	bdy = pos[1].y - pos[2].y;
	// TODO: check the sign here!
	vol = -0.5*(adx*bdy-bdx*ady);
	//vol = 0.5*(adx*bdy-bdx*ady);
#endif
	if(vol < 0.0) { // swap two vertices if the volume is negative
		swp = pos[0]; pos[0] = pos[1]; pos[1] = swp;
		swp = vel[0]; vel[0] = vel[1]; vel[1] = swp;
		vol *= -1;
	}
	return vol;
}


void psi_voxelize_annihilation(psi_rvec* pos0, psi_rvec* vel0, psi_real mass0, 
		psi_rvec* pos1, psi_rvec* vel1, psi_real mass1, psi_rvec* mbox, psi_grid* grid) {

	psi_int i, gridind;
	psi_real moments[10];
	psi_poly curpoly;
	psi_plane faces[PSI_NDIM+1];
	psi_voxels vox;

	// get the volume and correctly orient the tets
	// TODO: volume for degenerate tets?
	psi_real vol0 = psi_orient_tet(pos0, vel0);
	psi_real vol1 = psi_orient_tet(pos1, vel1);
	psi_real rho0 = mass0/vol0;
	psi_real rho1 = mass1/vol1;

	// initialize and clip to get the polyhedron
	psi_init_tet(&curpoly, pos0);
	psi_tet_faces_from_verts(faces, pos1);

	// voxelize the intersection
	psi_clip(&curpoly, faces, 4);
	psi_voxels_init(&vox, &curpoly, 0, mbox, grid);
	while(psi_voxels_next(&vox, moments, &gridind)) {
		grid->fields[PSI_GRID_M][gridind] += rho0*rho1*vol0*moments[0];
	}
}

#if 0

void psi_point_sample_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_dest_grid* grid) {

	psi_int i, gridind;
	psi_real cmin, cmax, dv, vol, bcoords[4];
	psi_rvec spos;
	psi_dvec ibox[2], ind;

	// clamp the index box to see if the poly needs to be clipped
	// against the window boundary 
	for(i = 0; i < PSI_NDIM; ++i) {
		cmin = rbox[0].xyz[i];
		cmax = rbox[1].xyz[i];
		if(grid->window[0][i] > cmin)
			cmin = grid->window[0][i];
		if(grid->periodic && grid->box[0][i] > cmin)
			cmin = grid->box[0][i];
		if(grid->window[1][i] < cmax)
			cmax = grid->window[1][i];
		if(grid->periodic && grid->box[1][i] < cmax)
			cmax = grid->box[1][i];
		ibox[0].ijk[i] = floor((cmin-grid->window[0][i])/grid->d[i]);
		ibox[1].ijk[i] = ceil((cmax-grid->window[0][i])/grid->d[i]);
		if(ibox[0].ijk[i] < 0) ibox[0].ijk[i] = 0;
		if(ibox[1].ijk[i] > grid->n[i]) ibox[1].ijk[i] = grid->n[i];
	}

#if PSI_NDIM == 3
	// voxel volume
	dv = grid->d[0]*grid->d[1]*grid->d[2];
	for(ind.i = ibox[0].i; ind.i < ibox[1].i; ++ind.i)
	for(ind.j = ibox[0].j; ind.j < ibox[1].j; ++ind.j)
	for(ind.k = ibox[0].k; ind.k < ibox[1].k; ++ind.k) {
	
		for(i = 0; i < PSI_NDIM; ++i)
			spos.xyz[i] = grid->window[0][i]+grid->d[i]*ind.ijk[i];

		// get the barycentric coordinates from the matrix inverse
		// return if the point lies outside any face of the tet
		vol = psi_barycentric(pos, spos, bcoords);
		if(!vol) goto next_sample;
		for(i = 0; i < PSI_NDIM+1; ++i) 
			if(bcoords[i] < 0.0) goto next_sample;

		gridind = grid->n[1]*grid->n[2]*ind.i + grid->n[2]*ind.j + ind.k;
		if(grid->fields[PSI_GRID_M]) grid->fields[PSI_GRID_M][gridind] += dv*mass/fabs(vol);

		// TODO: velocity stuff with BC coordinates 
	
		next_sample: continue;
	}
#elif PSI_NDIM == 2
	// voxel volume
	dv = grid->d[0]*grid->d[1];
	for(ind.i = ibox[0].i; ind.i < ibox[1].i; ++ind.i)
	for(ind.j = ibox[0].j; ind.j < ibox[1].j; ++ind.j) {
	
		for(i = 0; i < PSI_NDIM; ++i)
			spos.xyz[i] = grid->window[0][i]+grid->d[i]*ind.ijk[i];

		// get the barycentric coordinates from the matrix inverse
		// return if the point lies outside any face of the tet
		vol = psi_barycentric(pos, spos, bcoords);
		if(!vol) goto next_sample;
		for(i = 0; i < PSI_NDIM+1; ++i) 
			if(bcoords[i] < 0.0) goto next_sample;

		gridind = grid->n[1]*ind.i + ind.j;
		if(grid->fields[PSI_GRID_M]) grid->fields[PSI_GRID_M][gridind] += dv*mass/fabs(vol);

		// TODO: velocity stuff with BC coordinates 
	
		next_sample: continue;
	}

#endif
}

#endif

psi_real psi_barycentric(psi_rvec* pos, psi_rvec samppos, psi_real* bcoords) {

	psi_int i, j;
	psi_real btot;
	psi_real woff[PSI_NDIM+1];	
	psi_rvec fnorms[PSI_NDIM+1];	

	// Matrix inverse for the barycentric coordinates
#if PSI_NDIM == 3
	psi_real x1 = pos[0].x; psi_real y1 = pos[0].y; psi_real z1 = pos[0].z; 
	psi_real x2 = pos[1].x; psi_real y2 = pos[1].y; psi_real z2 = pos[1].z; 
	psi_real x3 = pos[2].x; psi_real y3 = pos[2].y; psi_real z3 = pos[2].z; 
	psi_real x4 = pos[3].x; psi_real y4 = pos[3].y; psi_real z4 = pos[3].z; 
	woff[0] = x2*(y3*z4 - y4*z3) + x3*(y4*z2 - y2*z4) + x4*(y2*z3 - y3*z2);
	woff[1] = x1*(y4*z3 - y3*z4) + x3*(y1*z4 - y4*z1) + x4*(y3*z1 - y1*z3);
	woff[2] = x1*(y2*z4 - y4*z2) + x2*(y4*z1 - y1*z4) + x4*(y1*z2 - y2*z1);
	woff[3] = x1*(y3*z2 - y2*z3) + x2*(y1*z3 - y3*z1) + x3*(y2*z1 - y1*z2);
	fnorms[0].x = ((y4 - y2)*(z3 - z2) - (y3 - y2)*(z4 - z2));
	fnorms[1].x = ((y3 - y1)*(z4 - z3) - (y3 - y4)*(z1 - z3));
	fnorms[2].x = ((y2 - y4)*(z1 - z4) - (y1 - y4)*(z2 - z4));
	fnorms[3].x = ((y1 - y3)*(z2 - z1) - (y1 - y2)*(z3 - z1));
	fnorms[0].y = ((x3 - x2)*(z4 - z2) - (x4 - x2)*(z3 - z2));
	fnorms[1].y = ((x4 - x3)*(z3 - z1) - (x1 - x3)*(z3 - z4));
	fnorms[2].y = ((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4));
	fnorms[3].y = ((x2 - x1)*(z1 - z3) - (x3 - x1)*(z1 - z2));
	fnorms[0].z = ((x4 - x2)*(y3 - y2) - (x3 - x2)*(y4 - y2));
	fnorms[1].z = ((x3 - x1)*(y4 - y3) - (x3 - x4)*(y1 - y3));
	fnorms[2].z = ((x2 - x4)*(y1 - y4) - (x1 - x4)*(y2 - y4));
	fnorms[3].z = ((x1 - x3)*(y2 - y1) - (x1 - x2)*(y3 - y1));
#elif PSI_NDIM == 2
	psi_real x1 = pos[0].x; psi_real y1 = pos[0].y; 
	psi_real x2 = pos[1].x; psi_real y2 = pos[1].y; 
	psi_real x3 = pos[2].x; psi_real y3 = pos[2].y; 
	woff[0] = x2*y3-x3*y2; 
	woff[1] = x3*y1-x1*y3;
	woff[2] = x1*y2-x2*y1;
	fnorms[0].x = y2-y3;
	fnorms[1].x = y3-y1;
	fnorms[2].x = y1-y2;
	fnorms[0].y = x3-x2;
	fnorms[1].y = x1-x3;
	fnorms[2].y = x2-x1;
#endif

	// renormalize such that bcoords sums to 1.0 
	btot = 0.0;
	for(i = 0; i < PSI_NDIM+1; ++i) {
		bcoords[i] = woff[i];
		for(j = 0; j < PSI_NDIM; ++j) bcoords[i] += samppos.xyz[j]*fnorms[i].xyz[j]; 
		btot += bcoords[i];	
	} 
	for(i = 0; i < PSI_NDIM+1; ++i) 
		bcoords[i] /= btot;
	return fabs(PSI_TET_FAC*btot); // return the volume. Might be zero, in which case the bcoords are NaN 
}

void psi_clip(psi_poly* poly, psi_plane* planes, psi_int nplanes) {

	// direct access to vertex buffer
	psi_vertex* vertbuffer = poly->verts; 
	psi_int* nverts = &poly->nverts; 
	if(*nverts <= 0) return;

	// variable declarations
	psi_int i, v, p, np, onv, vcur, vnext, vstart, 
			pnext, numunclipped;
	psi_real invw;

	// signed distances to the clipping plane
	psi_real sdists[POLYSZ];
	psi_real smin, smax;

	// for marking clipped vertices
	psi_int clipped[POLYSZ];

	// loop over each clip plane
	for(p = 0; p < nplanes; ++p) {

		// calculate signed distances to the clip plane
		onv = *nverts;
		smin = 1.0e30;
		smax = -1.0e30;
		memset(&clipped, 0, sizeof(clipped));
		for(v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d; 
			for(i = 0; i < PSI_NDIM; ++i)
				sdists[v] += vertbuffer[v].pos.xyz[i]*planes[p].n.xyz[i];
			// TODO: use integers to track whether to skip completely
			if(sdists[v] < smin) smin = sdists[v];
			if(sdists[v] > smax) smax = sdists[v];
			if(sdists[v] < 0.0) clipped[v] = 1;
		}

		// skip this face if the poly lies entirely on one side of it 
		if(smin >= 0.0) continue;
		if(smax <= 0.0) {
			*nverts = 0;
			return;
		}

		// check all edges and insert new vertices on the bisected edges 
		for(vcur = 0; vcur < onv; ++vcur) {
			if(clipped[vcur]) continue;
			for(np = 0; np < PSI_NDIM; ++np) {
				vnext = vertbuffer[vcur].pnbrs[np];
				if(!clipped[vnext]) continue;
				invw = 1.0/(sdists[vcur]-sdists[vnext]); // guaranteed to be nonzero
				for(i = 0; i < PSI_NDIM; ++i) {
					vertbuffer[*nverts].pos.xyz[i] = invw*(vertbuffer[vnext].pos.xyz[i]*sdists[vcur] - vertbuffer[vcur].pos.xyz[i]*sdists[vnext]);
					vertbuffer[*nverts].q.xyz[i] = invw*(vertbuffer[vnext].q.xyz[i]*sdists[vcur] - vertbuffer[vcur].q.xyz[i]*sdists[vnext]);
				}
#if PSI_NDIM == 3
				vertbuffer[*nverts].pnbrs[0] = vcur;
				vertbuffer[vcur].pnbrs[np] = *nverts;
#elif PSI_NDIM == 2
				vertbuffer[*nverts].pnbrs[1-np] = vcur;
				vertbuffer[*nverts].pnbrs[np] = -1;
				vertbuffer[vcur].pnbrs[np] = *nverts;
#endif
				(*nverts)++;
			}

		}

		// for each new vert, search around the faces for its new neighbors
		// and doubly-link everything
		for(vstart = onv; vstart < *nverts; ++vstart) {
			vcur = vstart;
			vnext = vertbuffer[vcur].pnbrs[0];
#if PSI_NDIM == 3
			do {
				for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
				vcur = vnext;
				pnext = (np+1)%3;
				vnext = vertbuffer[vcur].pnbrs[pnext];
			} while(vcur < onv);
			vertbuffer[vstart].pnbrs[2] = vcur;
			vertbuffer[vcur].pnbrs[1] = vstart;
#elif PSI_NDIM == 2
			if(vertbuffer[vstart].pnbrs[1] >= 0) continue;
			do {
				vcur = vnext;
				vnext = vertbuffer[vcur].pnbrs[0];
			} while(vcur < onv);
			vertbuffer[vstart].pnbrs[1] = vcur;
			vertbuffer[vcur].pnbrs[0] = vstart;
#endif
		}

		if(*nverts >= POLYSZ) 
			psi_printf("WARNING: Overflowed vertex buffer. nverts = %d >= %d\n", *nverts, POLYSZ);

		// go through and compress the vertex list, removing clipped verts
		// and re-indexing accordingly (reusing `clipped` to re-index everything)
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!clipped[v]) {
				vertbuffer[numunclipped] = vertbuffer[v];
				clipped[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) 
			for(np = 0; np < 3; ++np)
				vertbuffer[v].pnbrs[np] = clipped[vertbuffer[v].pnbrs[np]];
	}
}

void psi_tet_faces_from_verts(psi_plane* faces, psi_rvec* verts) {
#define dot(va, vb) (va.x*vb.x + va.y*vb.y + va.z*vb.z)
	// TODO: 2D vs 3D
	// NOTE: These vectors are not normalized! 
#if PSI_NDIM == 3
	psi_rvec tmpcent;
	faces[0].n.x = ((verts[3].y - verts[1].y)*(verts[2].z - verts[1].z) 
			- (verts[2].y - verts[1].y)*(verts[3].z - verts[1].z));
	faces[0].n.y = ((verts[2].x - verts[1].x)*(verts[3].z - verts[1].z) 
			- (verts[3].x - verts[1].x)*(verts[2].z - verts[1].z));
	faces[0].n.z = ((verts[3].x - verts[1].x)*(verts[2].y - verts[1].y) 
			- (verts[2].x - verts[1].x)*(verts[3].y - verts[1].y));
	//norm(faces[0].n);
	tmpcent.x = ONE_THIRD*(verts[1].x + verts[2].x + verts[3].x);
	tmpcent.y = ONE_THIRD*(verts[1].y + verts[2].y + verts[3].y);
	tmpcent.z = ONE_THIRD*(verts[1].z + verts[2].z + verts[3].z);
	faces[0].d = -dot(faces[0].n, tmpcent);

	faces[1].n.x = ((verts[2].y - verts[0].y)*(verts[3].z - verts[2].z) 
			- (verts[2].y - verts[3].y)*(verts[0].z - verts[2].z));
	faces[1].n.y = ((verts[3].x - verts[2].x)*(verts[2].z - verts[0].z) 
			- (verts[0].x - verts[2].x)*(verts[2].z - verts[3].z));
	faces[1].n.z = ((verts[2].x - verts[0].x)*(verts[3].y - verts[2].y) 
			- (verts[2].x - verts[3].x)*(verts[0].y - verts[2].y));
	//norm(faces[1].n);
	tmpcent.x = ONE_THIRD*(verts[2].x + verts[3].x + verts[0].x);
	tmpcent.y = ONE_THIRD*(verts[2].y + verts[3].y + verts[0].y);
	tmpcent.z = ONE_THIRD*(verts[2].z + verts[3].z + verts[0].z);
	faces[1].d = -dot(faces[1].n, tmpcent);

	faces[2].n.x = ((verts[1].y - verts[3].y)*(verts[0].z - verts[3].z) 
			- (verts[0].y - verts[3].y)*(verts[1].z - verts[3].z));
	faces[2].n.y = ((verts[0].x - verts[3].x)*(verts[1].z - verts[3].z) 
			- (verts[1].x - verts[3].x)*(verts[0].z - verts[3].z));
	faces[2].n.z = ((verts[1].x - verts[3].x)*(verts[0].y - verts[3].y) 
			- (verts[0].x - verts[3].x)*(verts[1].y - verts[3].y));
	//norm(faces[2].n);
	tmpcent.x = ONE_THIRD*(verts[3].x + verts[0].x + verts[1].x);
	tmpcent.y = ONE_THIRD*(verts[3].y + verts[0].y + verts[1].y);
	tmpcent.z = ONE_THIRD*(verts[3].z + verts[0].z + verts[1].z);
	faces[2].d = -dot(faces[2].n, tmpcent);

	faces[3].n.x = ((verts[0].y - verts[2].y)*(verts[1].z - verts[0].z) 
			- (verts[0].y - verts[1].y)*(verts[2].z - verts[0].z));
	faces[3].n.y = ((verts[1].x - verts[0].x)*(verts[0].z - verts[2].z) 
			- (verts[2].x - verts[0].x)*(verts[0].z - verts[1].z));
	faces[3].n.z = ((verts[0].x - verts[2].x)*(verts[1].y - verts[0].y) 
			- (verts[0].x - verts[1].x)*(verts[2].y - verts[0].y));
	//norm(faces[3].n);
	tmpcent.x = ONE_THIRD*(verts[0].x + verts[1].x + verts[2].x);
	tmpcent.y = ONE_THIRD*(verts[0].y + verts[1].y + verts[2].y);
	tmpcent.z = ONE_THIRD*(verts[0].z + verts[1].z + verts[2].z);
	faces[3].d = -dot(faces[3].n, tmpcent);
#elif PSI_NDIM == 2
#endif
}

void psi_split_coord(psi_poly* inpoly, psi_poly* outpolys, psi_real coord, psi_int ax) {

	psi_int v, np, i, npnxt, onv, vcur, vnext, vstart, pnext, nright;
	psi_real invw;
	psi_rvec newpos, newq;
	psi_int* nverts = &inpoly->nverts;
	psi_vertex* vertbuffer = inpoly->verts; 
	psi_int side[POLYSZ];
	psi_real sdists[POLYSZ];

	// calculate signed distances to the clip plane
	nright = 0;
	memset(&side, 0, sizeof(side));
	for(v = 0; v < *nverts; ++v) {
		sdists[v] = coord - vertbuffer[v].pos.xyz[ax];
		if(sdists[v] < 0.0) {
			side[v] = 1;
			nright++;
		}
	}

	// return if the poly lies entirely on one side of it 
	if(nright == 0) {
		outpolys[0] = *inpoly;
		outpolys[1].nverts = 0;
		return;
	}
	if(nright == *nverts) {
		outpolys[1] = *inpoly;
		outpolys[0].nverts = 0;
		return;
	}

	// check all edges and insert new vertices on the bisected edges 
	onv = inpoly->nverts;
	for(vcur = 0; vcur < onv; ++vcur) {
		if(side[vcur]) continue;
#if PSI_NDIM == 3
		for(np = 0; np < 3; ++np) {
			vnext = vertbuffer[vcur].pnbrs[np];
			if(!side[vnext]) continue;
			invw = 1.0/(sdists[vcur]-sdists[vnext]); // guaranteed to be nonzero
			for(i = 0; i < 3; ++i) {
				newpos.xyz[i] = invw*(vertbuffer[vnext].pos.xyz[i]*sdists[vcur] - vertbuffer[vcur].pos.xyz[i]*sdists[vnext]);
				newq.xyz[i] = invw*(vertbuffer[vnext].q.xyz[i]*sdists[vcur] - vertbuffer[vcur].q.xyz[i]*sdists[vnext]);
			}
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			vertbuffer[*nverts].pnbrs[0] = vcur;
			vertbuffer[vcur].pnbrs[np] = *nverts;
			(*nverts)++;
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			side[*nverts] = 1;
			vertbuffer[*nverts].pnbrs[0] = vnext;
			for(npnxt = 0; npnxt < 3; ++npnxt) 
				if(vertbuffer[vnext].pnbrs[npnxt] == vcur) break;
			vertbuffer[vnext].pnbrs[npnxt] = *nverts;
			(*nverts)++;
		}
#elif PSI_NDIM == 2
		for(np = 0; np < 2; ++np) {
			vnext = vertbuffer[vcur].pnbrs[np];
			if(!side[vnext]) continue;
			invw = 1.0/(sdists[vcur]-sdists[vnext]); // guaranteed to be nonzero
			for(i = 0; i < 2; ++i) {
				newpos.xyz[i] = invw*(vertbuffer[vnext].pos.xyz[i]*sdists[vcur] - vertbuffer[vcur].pos.xyz[i]*sdists[vnext]);
				newq.xyz[i] = invw*(vertbuffer[vnext].q.xyz[i]*sdists[vcur] - vertbuffer[vcur].q.xyz[i]*sdists[vnext]);
			}
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			vertbuffer[*nverts].pnbrs[1-np] = vcur;
			vertbuffer[*nverts].pnbrs[np] = -1;
			vertbuffer[vcur].pnbrs[np] = *nverts;
			(*nverts)++;
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			side[*nverts] = 1;
			vertbuffer[*nverts].pnbrs[np] = vnext;
			vertbuffer[*nverts].pnbrs[1-np] = -1;
			vertbuffer[vnext].pnbrs[1-np] = *nverts;
			(*nverts)++;
		}
#endif
	}

	// for each new vert, search around the faces for its new neighbors
	// and doubly-link everything
	for(vstart = onv; vstart < *nverts; ++vstart) {
		vcur = vstart;
		vnext = vertbuffer[vcur].pnbrs[0];
#if PSI_NDIM == 3
		do {
			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			vnext = vertbuffer[vcur].pnbrs[pnext];
		} while(vcur < onv);
		vertbuffer[vstart].pnbrs[2] = vcur;
		vertbuffer[vcur].pnbrs[1] = vstart;
#elif PSI_NDIM == 2
		if(vertbuffer[vstart].pnbrs[1] >= 0) continue;
		do {
			vcur = vnext;
			vnext = vertbuffer[vcur].pnbrs[0];
		} while(vcur < onv);
		vertbuffer[vstart].pnbrs[1] = vcur;
		vertbuffer[vcur].pnbrs[0] = vstart;
#endif
	}

	if(*nverts >= POLYSZ) 
		psi_printf("WARNING: Overflowed vertex buffer. nverts = %d >= %d\n", *nverts, POLYSZ);

	// copy and compress vertices into their new buffers
	// TODO: do this beforehand, since we have two buffers ready
	// reusing side[] for reindexing
	onv = *nverts;
	outpolys[0].nverts = 0;
	outpolys[1].nverts = 0;
	for(v = 0; v < onv; ++v) {
		outpolys[side[v]].verts[outpolys[side[v]].nverts] = vertbuffer[v];
		side[v] = outpolys[side[v]].nverts++;
	}
	for(v = 0; v < outpolys[0].nverts; ++v) 
		for(np = 0; np < PSI_NDIM; ++np)
			outpolys[0].verts[v].pnbrs[np] = side[outpolys[0].verts[v].pnbrs[np]];
	for(v = 0; v < outpolys[1].nverts; ++v) 
		for(np = 0; np < PSI_NDIM; ++np)
			outpolys[1].verts[v].pnbrs[np] = side[outpolys[1].verts[v].pnbrs[np]];
}
	
void psi_reduce(psi_poly* poly, psi_real* moments, psi_int polyorder, psi_int weight) {
	// reduce according to Powell and Abel (2015) for simplicity
#if PSI_NDIM == 3
	const static int nmom[3] = {1, 4, 10};
#elif PSI_NDIM == 2
	const static int nmom[3] = {1, 3, 6};
#endif
	psi_real mass;
	psi_int np, vstart, pstart, vcur, vnext, pnext;
	psi_rvec q0, q1, q2; 
	psi_vertex* vertbuffer = poly->verts; 
	psi_int* nverts = &poly->nverts; 
	psi_int emarks[*nverts][PSI_NDIM];
	memset(moments, 0, nmom[polyorder]*sizeof(psi_real)); 
	memset(&emarks, 0, sizeof(emarks));
#if PSI_NDIM == 3
	for(vstart = 0; vstart < *nverts; ++vstart)
	for(pstart = 0; pstart < 3; ++pstart) {
		if(emarks[vstart][pstart]) continue;
		pnext = pstart; 
		vcur = vstart; 
		emarks[vcur][pnext] = 1;
		vnext = vertbuffer[vcur].pnbrs[pnext];
		q0 = vertbuffer[vcur].q;
		for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
		vcur = vnext;
		pnext = (np+1)%3;
		emarks[vcur][pnext] = 1;
		vnext = vertbuffer[vcur].pnbrs[pnext];
		while(vnext != vstart) {
			q1 = vertbuffer[vnext].q;
			q2 = vertbuffer[vcur].q;

			mass = (-q2.x*q1.y*q0.z + q1.x*q2.y*q0.z + q2.x*q0.y*q1.z
			   	- q0.x*q2.y*q1.z - q1.x*q0.y*q2.z + q0.x*q1.y*q2.z); 
			moments[0] += mass;
			if(polyorder > 0) {
				moments[1] += 0.25*mass*(q0.x + q1.x + q2.x);
				moments[2] += 0.25*mass*(q0.y + q1.y + q2.y);
				moments[3] += 0.25*mass*(q0.z + q1.z + q2.z);
			}
			if(polyorder > 1) {
				moments[4] += 0.1*mass*(q0.x*q0.x + q1.x*q1.x + q2.x*q2.x + q1.x*q2.x + q0.x*(q1.x + q2.x));
				moments[5] += 0.05*mass*(q2.x*q0.y + q2.x*q1.y + 2*q2.x*q2.y + q0.x*(2*q0.y + q1.y + q2.y) + q1.x*(q0.y + 2*q1.y + q2.y));
				moments[6] += 0.05*mass*(q2.x*q0.z + q2.x*q1.z + 2*q2.x*q2.z + q0.x*(2*q0.z + q1.z + q2.z) + q1.x*(q0.z + 2*q1.z + q2.z));
				moments[7] += 0.1*mass*(q0.y*q0.y + q1.y*q1.y + q2.y*q2.y + q1.y*q2.y + q0.y*(q1.y + q2.y));
				moments[8] += 0.05*mass*(q2.y*q0.z + q2.y*q1.z + 2*q2.y*q2.z + q0.y*(2*q0.z + q1.z + q2.z) + q1.y*(q0.z + 2*q1.z + q2.z));
				moments[9] += 0.1*mass*(q0.z*q0.z + q1.z*q1.z + q2.z*q2.z + q1.z*q2.z + q0.z*(q1.z + q2.z));
			}
			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			emarks[vcur][pnext] = 1;
			vnext = vertbuffer[vcur].pnbrs[pnext];
		}
	}	
#elif PSI_NDIM == 2
	for(vcur = 0; vcur < *nverts; ++vcur) {
		vnext = vertbuffer[vcur].pnbrs[0];
		q0 = vertbuffer[vcur].q;
		q1 = vertbuffer[vnext].q;
		mass = (q0.x*q1.y - q0.y*q1.x); 
		moments[0] += mass;
		if(polyorder > 0) {
			moments[1] += ONE_THIRD*mass*(q0.x + q1.x);
			moments[2] += ONE_THIRD*mass*(q0.y + q1.y);
		}
		if(polyorder > 1) {
			moments[3] += ONE_SIXTH*mass*(q0.x*q0.x + q1.x*q1.x + q0.x*q1.x);
			moments[4] += ONE_TWLTH*mass*(q0.x*(2*q0.y + q1.y) + q1.x*(q0.y + 2*q1.y));
			moments[5] += ONE_SIXTH*mass*(q0.y*q0.y + q1.y*q1.y + q0.y*q1.y);
		}
	}	
#endif

}

void psi_init_tet(psi_poly* poly, psi_rvec* pos) {
	// initialize the tet as an edge-vertex graph 
	// fill in mass coordinates too
	psi_int v;
	memset(poly, 0, sizeof(psi_poly));
#if PSI_NDIM == 3
	poly->nverts = 4;
	poly->verts[0].pnbrs[0] = 1;	
	poly->verts[0].pnbrs[1] = 3;	
	poly->verts[0].pnbrs[2] = 2;	
	poly->verts[1].pnbrs[0] = 2;	
	poly->verts[1].pnbrs[1] = 3;	
	poly->verts[1].pnbrs[2] = 0;	
	poly->verts[2].pnbrs[0] = 0;	
	poly->verts[2].pnbrs[1] = 3;	
	poly->verts[2].pnbrs[2] = 1;	
	poly->verts[3].pnbrs[0] = 1;	
	poly->verts[3].pnbrs[1] = 2;	
	poly->verts[3].pnbrs[2] = 0;	
	// TODO: all four BC coordinates for symmetry?
	poly->verts[1].q.x = 1.0;
	poly->verts[2].q.y = 1.0;
	poly->verts[3].q.z = 1.0;
#elif PSI_NDIM == 2
	poly->nverts = 3;
	poly->verts[0].pnbrs[0] = 1;	
	poly->verts[0].pnbrs[1] = 2;	
	poly->verts[1].pnbrs[0] = 2;	
	poly->verts[1].pnbrs[1] = 0;	
	poly->verts[2].pnbrs[0] = 0;	
	poly->verts[2].pnbrs[1] = 1;	
	poly->verts[1].q.x = 1.0;
	poly->verts[2].q.y = 1.0;
#endif
	for(v = 0; v < PSI_NDIM+1; ++v) poly->verts[v].pos = pos[v];
}

void psi_voxels_init(psi_voxels* vox, psi_poly* poly, psi_int polyorder, psi_rvec* rbox, psi_grid* grid) {

	psi_int i, cliplo, cliphi, nplanes;
	psi_real cmin, cmax;
	psi_plane planes[2*PSI_NDIM];

	// first see if the poly needs to be clipped against the window boundary 
	// and compute the index bounds for the polyhedron
	memset(planes, 0, sizeof(planes));
	nplanes = 0;
	for(i = 0; i < PSI_NDIM; ++i) {
		cliplo = 0; cliphi = 0;
		cmin = rbox[0].xyz[i];
		cmax = rbox[1].xyz[i];
		if(grid->window[0].xyz[i] > cmin) {
			cmin = grid->window[0].xyz[i];
			cliplo = 1;
		}	
		if(grid->window[1].xyz[i] < cmax) {
			cmax = grid->window[1].xyz[i];
			cliphi = 1;
		}	
		if(cliplo) {
			planes[nplanes].n.xyz[i] = 1.0;
			planes[nplanes].d = -cmin;
			nplanes++;
		}
		if(cliphi) {
			planes[nplanes].n.xyz[i] = -1.0;
			planes[nplanes].d = cmax;
			nplanes++;
		}
		poly->ibox[0].ijk[i] = floor((cmin-grid->window[0].xyz[i])/grid->d.xyz[i]);
		poly->ibox[1].ijk[i] = ceil((cmax-grid->window[0].xyz[i])/grid->d.xyz[i]);
		if(poly->ibox[0].ijk[i] < 0) poly->ibox[0].ijk[i] = 0;
		if(poly->ibox[1].ijk[i] > grid->n.ijk[i]) poly->ibox[1].ijk[i] = grid->n.ijk[i];
	}
	psi_clip(poly, planes, nplanes);

	// initialize the stack
	vox->polyorder = polyorder;
	vox->grid = grid;
	vox->stack[0] = *poly;
	vox->nstack = 1;
}

psi_int psi_voxels_next(psi_voxels* vox, psi_real* moments, psi_int* gridind) {

	psi_int i, dmax, spax, siz;
	psi_poly curpoly;

	// pointers for easier acces
	psi_poly* stack = vox->stack;
	psi_int* nstack = &vox->nstack;
	psi_grid* grid = vox->grid;

	while(*nstack > 0) {
		curpoly = stack[--*nstack];
		if(curpoly.nverts <= 0) continue;
		dmax = 0; spax = 0;
		for(i = 0; i < PSI_NDIM; ++i) {
			siz = curpoly.ibox[1].ijk[i]-curpoly.ibox[0].ijk[i];
			if(siz > dmax) {
				dmax = siz; 
				spax = i;
			}	
		}
		if(dmax == 1) {
			// if all three axes are only one voxel long, reduce the single voxel to the dest grid
			// reduce and add it to the target grid if its mass is nonzero
			psi_reduce(&curpoly, moments, vox->polyorder, 0);
			if(moments[0] <= 0.0) continue;
#if PSI_NDIM == 3
			*gridind = grid->n.ijk[1]*grid->n.ijk[2]*curpoly.ibox[0].i + grid->n.ijk[2]*curpoly.ibox[0].j + curpoly.ibox[0].k;
#elif PSI_NDIM == 2
			*gridind = grid->n.ijk[1]*curpoly.ibox[0].i + curpoly.ibox[0].j;
#endif
			return 1;
		}

		// split the poly and push children to the stack
		psi_split_coord(&curpoly, &stack[*nstack], 
				grid->d.xyz[spax]*(curpoly.ibox[0].ijk[spax]+dmax/2)+grid->window[0].xyz[spax], spax);
		memcpy(stack[*nstack].ibox, curpoly.ibox, 2*sizeof(psi_dvec));
		stack[*nstack].ibox[1].ijk[spax] -= dmax-dmax/2; 
		memcpy(stack[*nstack+1].ibox, curpoly.ibox, 2*sizeof(psi_dvec));
		stack[*nstack+1].ibox[0].ijk[spax] += dmax/2;
		*nstack += 2;
	}
	return 0;
}

// gives the solid angle of the triangular cone
psi_real psi_omega3(psi_rvec v1, psi_rvec v2, psi_rvec v3) {
	psi_int i;
	psi_real det, div;

	// normal.ze v1, v2, v3
	div = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
	for(i = 0; i < 3; ++i) v1.xyz[i] /= div;
	div = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
	for(i = 0; i < 3; ++i) v2.xyz[i] /= div;
	div = sqrt(v3.x*v3.x + v3.y*v3.y + v3.z*v3.z);
	for(i = 0; i < 3; ++i) v3.xyz[i] /= div;

	// Solid angle of a triangle .y Oosterom and Strackee
	det = v1.x*(v2.y*v3.z-v3.y*v2.z)-v1.y*(v2.x*v3.z-v3.x*v2.z)+v1.z*(v2.x*v3.y-v3.x*v2.y);
	div = 1.0;
	for(i = 0; i < 3; ++i) {
		div += v1.xyz[i]*v2.xyz[i];
		div += v2.xyz[i]*v3.xyz[i];
		div += v3.xyz[i]*v1.xyz[i];
	}
	return 2.0*atan2(det, div);
} 



#if 0
psi_real psi_jacobdet_from_3cube(psi_rvec* in, psi_dvec ind) {

	// does the full triquadratic interpolation for now
	const static psi_real diff[3][3] = {{-3.0, 4.0, -1.0}, {-1.0, 0.0, 1.0}, {1.0, -4.0, 3.0}};
	psi_int i, j, k, ax;
	psi_real jacob[3][3];
	memset(jacob, 0, sizeof(jacob));
	for(ax = 0; ax < PSI_NDIM; ++ax) {
#if PSI_NDIM == 3
		for(i = 0; i < 3; ++i)
			jacob[ax][0] += diff[i][0]*in[ind3(i, ind.j, ind.k)].xyz[ax];
		for(j = 0; j < 3; ++j)
			jacob[ax][1] += diff[j][1]*in[ind3(ind.i, j, ind.k)].xyz[ax];
		for(k = 0; k < 3; ++k)
			jacob[ax][2] += diff[k][2]*in[ind3(ind.i, ind.j, k)].xyz[ax];
#elif PSI_NDIM == 2
#endif
	}
	return jacob[0][0]*((jacob[1][1]*jacob[2][2]) - (jacob[2][1]*jacob[1][2])) 
			-jacob[0][1]*(jacob[1][0]*jacob[2][2] - jacob[2][0]*jacob[1][2]) 
			+ jacob[0][2]*(jacob[1][0]*jacob[2][1] - jacob[2][0]*jacob[1][1]);

}
#endif

