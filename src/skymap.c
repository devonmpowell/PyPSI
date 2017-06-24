#include "skymap.h"


// the main function
// traces beams as defined by the grb_params passed in
void psi_skymap(psi_grid* grid, psi_mesh* mesh, psi_int bstep) {

	psi_int p, b, ax, e, v, tind, t; 
	psi_dvec grind;
	psi_real vol;
	psi_real moments[10];
	psi_rvec corners[64], center, brbox[2];
	psi_rvec gpos[mesh->dim+1], grbox[2];
	psi_int vpere = mesh->elemtype;
	psi_beam beam;
	psi_plane beamfaces[16];
	psi_tet_buffer tetbuf;

	psi_rvec rawverts[3]; 
	psi_rvec beamverts[6];

	psi_rvec tpos[vpere], tvel[vpere], trbox[2];
	psi_real tmass;

	if(grid->type != PSI_GRID_HPRING)
		return;

	psi_rvec obspos = {{20.,20.,20.}};

	psi_real mtot = 0.0;

	// First, build an rtree to get logN spatial queries
	psi_rtree rtree;
	psi_rtree_query qry;	
	psi_rtree_init(&rtree, 2*mesh->nelem);
	for(e = 0; e < mesh->nelem; ++e) {
		// make it periodic, get its bounding box, and check it against the grid
		for(v = 0; v < vpere; ++v) {
			tpos[v] = mesh->pos[mesh->connectivity[e*vpere+v]];
			mtot += (1.0/vpere)*mesh->mass[mesh->connectivity[e*vpere+v]];
		}

		// TODO: mesh->box may not exist for nonperiodic meshes!!!
		if(!psi_aabb_periodic(tpos, trbox, mesh->box, mesh)) continue;
		for(ax = 0; ax < 3; ++ax) {
			trbox[0].xyz[ax] -= obspos.xyz[ax];
			trbox[1].xyz[ax] -= obspos.xyz[ax];
		}
		// insert each element into the tree
		psi_rtree_insert(&rtree, trbox, e);
	}
	//psi_rtree_print(stdout, &rtree);
	psi_tet_buffer_init(&tetbuf, 1000.0, 0);

	printf("Total mass loaded = %.5e\n", mtot);

				/// TODO: test
				mtot = 0.0;
				for(ax = 0; ax < 3; ++ax) {
					brbox[0].xyz[ax] = -40.0;
					brbox[1].xyz[ax] = 40.0;
				}
				psi_rtree_query_init(&qry, &rtree, brbox);
				while(psi_rtree_query_next(&qry, &e)) {

					tmass = 0.0;
					for(v = 0; v < vpere; ++v) {
						tind = mesh->connectivity[e*vpere+v];
						tpos[v] = mesh->pos[tind];
						for(ax = 0; ax < 3; ++ax)
							tpos[v].xyz[ax] -= obspos.xyz[ax];
						tvel[v] = mesh->vel[tind];
						tmass += (1.0/vpere)*mesh->mass[tind]; // TODO: hack
					}

					mtot += tmass;




				}
	printf("Total mass queried = %.5e\n", mtot);

	

	// for every pixel in the image...
	memset(grid->fields[0], 0, grid->n.k*sizeof(psi_real));
	for(p = 0; p < grid->n.k; ++p) {

		// get the four corner vectors and the center
		// this depends on the pixelization scheme given by grbp
		grind.i = p;
		if(!psi_grid_get_cell_geometry(grid, grind, bstep, corners, &center, &vol)) 
			continue;

		// divide the pixel into four triangular beams,
		// all four share the center point 
		for(b = 0; b < 4*bstep; ++b) {
	
			// set spatial indices of beam four-vectors, starting at index 1
			for(ax = 0; ax < 3; ++ax) {
				rawverts[0].xyz[ax] = corners[b].xyz[ax];
				rawverts[1].xyz[ax] = corners[(b+1)%(4*bstep)].xyz[ax];
				rawverts[2].xyz[ax] = center.xyz[ax]; 
			}
	
			// TODO: rotate/transform/boost the beam position to the source
	
			// trace the beam, then accumulate the
			// resulting flux into the pixel
			//trace_beam(grbp, &beam);

			psi_real rmin = 0.01;
			psi_real rmax = 200.0;
			psi_real rstep = 0.5;
			psi_real rold = rmin;
			psi_real rnew;
			while(rold < rmax) {
				rnew = rold+rstep;

				// find the bounding box and query all elements that
				// may intersect it
				for(v = 0; v < 3; ++v) 
				for(ax = 0; ax < 3; ++ax) {
					beamverts[v+0].xyz[ax] = rold*rawverts[v].xyz[ax];
					beamverts[v+3].xyz[ax] = rnew*rawverts[v].xyz[ax];
				} 
				psi_aabb(beamverts, 6, brbox);

				/// TODO: test
				//for(ax = 0; ax < 3; ++ax) {
					//brbox[0].xyz[ax] -= 10.0;
					//brbox[1].xyz[ax] += 10.0;
				//}

				// get the beam faces
				for(v = 0; v < 3; ++v) {
					psi_cross3(beamfaces[v].n, rawverts[v], rawverts[(v+1)%3]);
					beamfaces[v].d = 0.0; 
				}
				psi_rvec tmp0, tmp1;
				for(ax = 0; ax < 3; ++ax) {
					tmp0.xyz[ax] = beamverts[0].xyz[ax] - beamverts[2].xyz[ax];
					tmp1.xyz[ax] = beamverts[1].xyz[ax] - beamverts[2].xyz[ax];
				}
				psi_cross3(beamfaces[3].n, tmp0, tmp1); 
				beamfaces[3].d = -psi_dot3(beamfaces[3].n, beamverts[2]);
				for(ax = 0; ax < 3; ++ax) {
					tmp0.xyz[ax] = beamverts[3].xyz[ax] - beamverts[5].xyz[ax];
					tmp1.xyz[ax] = beamverts[4].xyz[ax] - beamverts[5].xyz[ax];
				}
				psi_cross3(beamfaces[4].n, tmp1, tmp0); 
				beamfaces[4].d = -psi_dot3(beamfaces[4].n, beamverts[5]);


				psi_rtree_query_init(&qry, &rtree, brbox);
				while(psi_rtree_query_next(&qry, &e)) {

					// a local copy of the element, in case it's modified
					// make it periodic, get its bounding box, and check it against the grid
					tmass = 0.0;
					for(v = 0; v < vpere; ++v) {
						tind = mesh->connectivity[e*vpere+v];
						tpos[v] = mesh->pos[tind];
						for(ax = 0; ax < 3; ++ax)
							tpos[v].xyz[ax] -= obspos.xyz[ax];
						tvel[v] = mesh->vel[tind];
						tmass += (1.0/vpere)*mesh->mass[tind]; // TODO: hack
					}

					// refine elements into the tet buffer 
					// loop over each tet in the buffer
					psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, mesh->elemtype);
					for(t = 0; t < tetbuf.num; ++t) {
			
						// copy the position to the ghost array and compute its aabb
						memcpy(gpos, &tetbuf.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
						//psi_aabb(gpos, mesh->dim+1,grbox);
						// TODO: check the boxes before clipping....
						psi_clip_reduce_tet(gpos, beamfaces, 5, moments);
						if(moments[0] <= 0.0) continue;

						grid->fields[0][p] += tetbuf.mass[t]*moments[0]; 

						//printf("clip and reduce; mass = %.5e\n", moments[0]);
			
						// make ghosts and sample tets to the grid
						//psi_make_ghosts(gpos, grbox, &nghosts, (mesh->dim+1), grid->window, mesh);
						//for(g = 0; g < nghosts; ++g) {
							//psi_rvec* curtet = &gpos[(mesh->dim+1)*g];
							//if(grid->sampling == PSI_SAMPLING_VOLUME)
								//psi_voxelize_tet(gpos, &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], grbox);
							//else if(grid->sampling == PSI_SAMPLING_POINT)
								//psi_point_sample_tet(&gpos[(mesh->dim+1)*g], &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
						//}
					}

				
				}
			
				rold = rnew;
			}

		}
	}

	psi_rtree_destroy(&rtree);

#if 0
		// refine elements into the tet buffer 
		// loop over each tet in the buffer
		psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, mesh->elemtype);
		for(t = 0; t < tetbuf.num; ++t) {

			// copy the position to the ghost array and compute its aabb
			memcpy(gpos, &tetbuf.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
			psi_aabb(gpos, mesh->dim+1,grbox);

			// make ghosts and sample tets to the grid
			psi_make_ghosts(gpos, grbox, &nghosts, (mesh->dim+1), grid->window, mesh);
			for(g = 0; g < nghosts; ++g) {
				psi_rvec* curtet = &gpos[(mesh->dim+1)*g];
				//if(grid->sampling == PSI_SAMPLING_VOLUME)
					psi_voxelize_tet(&gpos[(mesh->dim+1)*g], &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
				//else if(grid->sampling == PSI_SAMPLING_POINT)
					//psi_point_sample_tet(&gpos[(mesh->dim+1)*g], &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
			}
		}
#endif


	// if we are in healpix, solve for the correction factors
	// that make it truly equal-area
	//if(grbp->pix_type == PIX_HEALPIX)
		//correct_healpix_edges(grbp, 1.0e-6);
} 

#if 0
// traces a beam defined by the current 
void trace_beam(grb_params* grbp, grb_beam* beam) {

	int v, ax;
	real r2;
	grb_beam btmp;

	// make a local copy
	btmp = *beam;
	btmp.flux = 0.0;

	if(grbp->pix_type == PIX_ORTHO) {
		int andcmp = 1;
		int orcmp = 0;
		for(v = 0; v < 3; ++v) {
		
			r2 = btmp.pos[v][1]*btmp.pos[v][1]+btmp.pos[v][2]*btmp.pos[v][2];
			andcmp &= (r2 < 1.0);
			orcmp |= (r2 < 1.0);
		}
	
		if(andcmp)
			btmp.flux = 1.0;
	
	}
	else if(grbp->pix_type == PIX_HEALPIX) {
		//btmp.flux = btmp.pos[2][3];
		//
		// spatial indices start at 1!
		btmp.flux = omega(&btmp.pos[0][1], &btmp.pos[1][1], &btmp.pos[2][1]); 
	}



	*beam = btmp;
}

void correct_healpix_edges(grb_params* grbp, real tol) {

		printf("Healpix detected. Solving the correction matrix for curved edges...\n");

		int p, n, b, ax, iter;
		int neighb[8];
		real max_err, err;
		real corners[4][3];
		real center[3];
		real ptmp[3][3];
		real* transfer, *allomega;

		// allocate an array for the (fractional) flux transfers between edges
		// identified by a pixel and its eastern neighbors
		transfer = malloc(2*grbp->npix*sizeof(real));
		allomega = malloc(grbp->npix*sizeof(real));
		memset(transfer, 0, 2*grbp->npix*sizeof(real));
		memset(allomega, 0, grbp->npix*sizeof(real));

		// get the solid angle as quoted by flat pixel edges 
		// divide the pixel into four triangular beams,
		// all four share the center point 
		for(p = 0; p < grbp->npix; ++p) {
			get_pixel_vectors(grbp, p, &corners[0][0], center);
			for(b = 0; b < 4; ++b) {
				for(ax = 0; ax < 3; ++ax) {
					ptmp[0][ax] = corners[b][ax];
					ptmp[1][ax] = corners[(b+1)%4][ax];
					ptmp[2][ax] = center[ax]; 
				}
				allomega[p] += omega(ptmp[0], ptmp[1], ptmp[2])*grbp->npix/FOUR_PI; 
			}
		}


		// iteratively relax the grid to find the solution
		max_err = 1.0e10;
		for(iter = 0; iter < 10 && max_err > tol; ++iter) {

			max_err = 0.0;
			for(p = 0; p < grbp->npix; ++p) {
	
				// get the error
				err = fabs(1.0-allomega[p]);
				if(err > max_err) max_err = err;
	
				// get the pixel neighbors
				neighbors(grbp, p, neighb);
				printf("Got neighbors for pixel %d:\n", p);
				for(n = 0; n < 8; ++n)
					printf(" - %d\n", neighb[n]);
				printf(" - total solid angle = %.5e, err = %.5e\n", allomega[p], err); 
			
			}
		
			printf(" -- Iter %d, max_err = %.5e\n", iter, max_err);
		
		}


		free(transfer);
		free(allomega);

}

#endif
