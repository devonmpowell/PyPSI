#include "skymap.h"


// the main function
// traces beams as defined by the grb_params passed in
void psi_skymap(psi_grid* grid, psi_mesh* mesh, psi_int bstep, psi_int mode) {

	psi_int p, b, ax, v, tind, t, e0, e1, t0, t1, m; 
	psi_dvec grind;
	psi_real vol;
	psi_real moments[10];
	psi_real vertweights[4];
	psi_rvec corners[64], center, brbox[2];
	psi_rvec gpos[mesh->dim+1], grbox[2], gpos0[mesh->dim+1], gpos1[mesh->dim+1], grbox0[2], grbox1[2], gvel0[mesh->dim+1], gvel1[mesh->dim+1];
	psi_int vpere = mesh->elemtype;
	psi_plane beamfaces[16];
	psi_tet_buffer tetbuf0, tetbuf1;

								psi_poly curpoly, tpoly;
									psi_plane tetfaces[4];

	psi_rvec rawverts[3]; 
	psi_rvec beamverts[6];

	psi_rvec tpos0[vpere], tvel0[vpere], trbox0[2], erbox0[2], erbox1[2], qbox0[2], qbox1[2], qbox2[2], qbox3[2];
	psi_rvec tpos1[vpere], tvel1[vpere], trbox1[2];
	psi_real tmass0, tmass1;

	if(grid->type != PSI_GRID_HPRING)
		return;

	psi_int nfail = 0;

	psi_rvec obspos = {{20.,20.,20.}};
	//psi_rvec obspos = {{0.,0.,0.}};

	psi_real mtot = 0.0;

	// First, build an rtree to get logN spatial queries
	psi_rtree rtree;
	psi_rtree_query qry0, qry1;	
	psi_rtree_init(&rtree, 2*mesh->nelem);
	for(e0 = 0; e0 < mesh->nelem; ++e0) {
		// make it periodic, get its bounding box, and check it against the grid
		for(v = 0; v < vpere; ++v) {
			tpos0[v] = mesh->pos[mesh->connectivity[e0*vpere+v]];
			mtot += (1.0/vpere)*mesh->mass[mesh->connectivity[e0*vpere+v]];
		}

		// TODO: mesh->box may not exist for nonperiodic meshes!!!
		if(!psi_aabb_periodic(tpos0, trbox0, mesh->box, mesh)) continue;
		for(ax = 0; ax < 3; ++ax) {
			trbox0[0].xyz[ax] -= obspos.xyz[ax];
			trbox0[1].xyz[ax] -= obspos.xyz[ax];
		}
		// insert each element into the tree
		psi_rtree_insert(&rtree, trbox0, e0);
	}
	psi_tet_buffer_init(&tetbuf0, 1000.0, 0);
	psi_tet_buffer_init(&tetbuf1, 1000.0, 0);
	//psi_tet_buffer_init(&tetbuf0, 1.0, 4);
	//psi_tet_buffer_init(&tetbuf1, 1.0, 4);

	/// RTREE test
	// see whether all bound boxes completely contain their children
	//psi_rtree_print(stdout, &rtree);
	printf("Total mass loaded = %.5e\n", mtot);
	/// TODO: test
	mtot = 0.0;
	for(ax = 0; ax < 3; ++ax) {
		brbox[0].xyz[ax] = mesh->box[0].xyz[ax] - obspos.xyz[ax]; 
		brbox[1].xyz[ax] = mesh->box[1].xyz[ax] - obspos.xyz[ax]; 
	}
	psi_rtree_query_init(&qry0, &rtree, brbox);
	while(psi_rtree_query_next(&qry0, &e0)) {
		tmass0 = 0.0;
		for(v = 0; v < vpere; ++v) {
			tind = mesh->connectivity[e0*vpere+v];
			tmass0 += (1.0/vpere)*mesh->mass[tind]; 
		}
		mtot += tmass0;
	}
	printf("Total mass queried = %.5e\n", mtot);

	psi_real omtot = 0.0;
	psi_int* elemmarks = psi_malloc(mesh->nelem*sizeof(psi_int));
	memset(elemmarks, 0, mesh->nelem*sizeof(psi_int));
	mtot = 0.0;

	// for every pixel in the image...
	memset(grid->fields[0], 0, grid->n.k*sizeof(psi_real));
	for(p = 0; p < grid->n.k; ++p) {

		// get the four corner vectors and the center
		// this depends on the pixelization scheme given by grbp
		grind.i = p;
		if(!psi_grid_get_cell_geometry(grid, grind, bstep, corners, &center, &vol)) {
			printf("Bad cell geometry!\n");
			continue;
		} 

		// divide the pixel into four triangular beams,
		// all four share the center point 
		for(b = 0; b < 4*bstep; ++b) {
	
			// set spatial indices of beam four-vectors, starting at index 1
			for(ax = 0; ax < 3; ++ax) {
				rawverts[0].xyz[ax] = corners[b].xyz[ax];
				rawverts[1].xyz[ax] = corners[(b+1)%(4*bstep)].xyz[ax];
				rawverts[2].xyz[ax] = center.xyz[ax]; 
			}
	
			//psi_real rmin = 10.0;
			//psi_real rmax = 10.0 + 0.001;
			//psi_real rstep = .001001; 
			

			psi_real rmin = 0.0; 
			psi_real rmax = 100.0; 
			psi_real rstep = 1.0; 
			
			psi_real rold = rmin;
			psi_real rnew;
			while(rold < rmax) {

				// find the bounding box and query all elements that
				// may intersect it
				rnew = rold+rstep;
				for(v = 0; v < 3; ++v) {
					for(ax = 0; ax < 3; ++ax) {
						beamverts[v+0].xyz[ax] = rold*rawverts[v].xyz[ax];
						beamverts[v+3].xyz[ax] = rnew*rawverts[v].xyz[ax];
					} 
				} 
				psi_aabb(beamverts, 6, brbox);

				// create the beam faces
				// TODO: can we do this outside of the loop???
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

				// use the r*-tree to query all elements that might
				// intersect the beam segment's aabb
				psi_rtree_query_init(&qry0, &rtree, brbox);
				while(psi_rtree_query_next(&qry0, &e0)) {

					// load the element pos, vel, mass 
					// refine elements into a buffer of linear tets 
					elemmarks[e0] = 1;
					tmass0 = 0.0;
					for(v = 0; v < vpere; ++v) {
						tind = mesh->connectivity[e0*vpere+v];
						tpos0[v] = mesh->pos[tind];
						for(ax = 0; ax < 3; ++ax)
							tpos0[v].xyz[ax] -= obspos.xyz[ax];
						tvel0[v] = mesh->vel[tind];
						tmass0 += (1.0/vpere)*mesh->mass[tind]; // TODO: hack
					}
					psi_tet_buffer_refine(&tetbuf0, tpos0, tvel0, tmass0, mesh->elemtype);

					// linear skymaps
					// signal proportional to the density
					if(mode == PSI_SKYMAP_RHO_LINEAR) {

						// loop over all linear tets in the refinement buffer
						for(t = 0; t < tetbuf0.num; ++t) {
				
							// copy the position to the ghost array and compute its aabb
							// skip if the aabb does not intersect the beam
							memcpy(gpos, &tetbuf0.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
							psi_aabb(gpos, mesh->dim+1, grbox);
							if(!psi_aabb_ixn(brbox, grbox, qbox2))
								continue;

							// clip and reduce the tet against the beam
							// skip if there is zero mass left
							psi_init_tet(&curpoly, gpos);
							psi_clip(&curpoly, beamfaces, 5);
							psi_reduce(&curpoly, moments, 1, 0);
							if(moments[0] <= 0.0) continue;
	
							// set up the vertex weights,
							// depending on the specified weight scheme
							for(v = 0; v < 4; ++v) {

								// mass-weighted for debugging
								vertweights[v] = tetbuf0.mass[t]; 
	
								// inverse r^2, gpos is already relative to the observer position
								//vertweights[v] = tetbuf0.mass[t]/(gpos[v].x*gpos[v].x+gpos[v].y*gpos[v].y+gpos[v].z*gpos[v].z); 
							}
	
							// subtract. Now moments[0->3] contain barycentric moments,
							// properly weighted wrt the original tet vertices
							// add the mass to the grid, linearly interpolating the vertex weights 
							moments[0] -= (moments[1]+moments[2]+moments[3]); 
							for(v = 0; v < 4; ++v) {
								grid->fields[0][p] += vertweights[v]*moments[v]; 
								mtot += vertweights[v]*moments[v]; 
							}					
						}
					}

					// quadratic skymaps
					// signal proportional to the squared density
					else if(mode == PSI_SKYMAP_RHO_SQUARED) {

						// get the element's aabb, find its intersection with the beams's aabb
						// then make a new query box out of the intersection
						psi_aabb(tpos0, vpere, erbox0);
						if(!psi_aabb_ixn(brbox, erbox0, qbox0))
							continue;

						// query elements close to both the beam and e0
						psi_rtree_query_init(&qry1, &rtree, qbox0);
						while(psi_rtree_query_next(&qry1, &e1)) {

							tmass1 = 0.0;
							for(v = 0; v < vpere; ++v) {
								tind = mesh->connectivity[e1*vpere+v];
								tpos1[v] = mesh->pos[tind];
								for(ax = 0; ax < 3; ++ax)
									tpos1[v].xyz[ax] -= obspos.xyz[ax];
								tvel1[v] = mesh->vel[tind];
								tmass1 += (1.0/vpere)*mesh->mass[tind]; // TODO: hack
							}

							// get the element's aabb, find its intersection with the beams's aabb
							// then make a new query box out of the intersection
							psi_aabb(tpos1, vpere, erbox1);
							if(!psi_aabb_ixn(qbox0, erbox1, qbox1)) 
								continue;

							// refine the second element now
							// and loop over pairs of linear tets (this bit is n^2)
							psi_tet_buffer_refine(&tetbuf1, tpos1, tvel1, tmass1, mesh->elemtype);
							for(t0 = 0; t0 < tetbuf0.num; ++t0) {
	
								// copy the position to the ghost array and compute its aabb
								memcpy(gpos0, &tetbuf0.pos[(mesh->dim+1)*t0], (mesh->dim+1)*sizeof(psi_rvec));
								memcpy(gvel0, &tetbuf0.vel[(mesh->dim+1)*t0], (mesh->dim+1)*sizeof(psi_rvec));
								psi_real mass0 = tetbuf0.mass[t0];
								psi_aabb(gpos0, mesh->dim+1, grbox0);
								if(!psi_aabb_ixn(qbox1, grbox0, qbox2))
									continue;

								// initialize tet0 and clip it against the beam segment
								psi_init_tet(&curpoly, gpos0);
								psi_clip(&curpoly, beamfaces, 5);
								if(curpoly.nverts <= 0) continue;
								tpoly = curpoly;
						
								// loop over all tets in buffer 1
								for(t1 = 0; t1 < tetbuf1.num; ++t1) {

									// copy the position to the ghost array and compute its aabb
									memcpy(gpos1, &tetbuf1.pos[(mesh->dim+1)*t1], (mesh->dim+1)*sizeof(psi_rvec));
									memcpy(gvel1, &tetbuf1.vel[(mesh->dim+1)*t1], (mesh->dim+1)*sizeof(psi_rvec));
									psi_real mass1 = tetbuf1.mass[t1];
									psi_aabb(gpos1, mesh->dim+1, grbox1);
									if(!psi_aabb_ixn(qbox2, grbox1, qbox3))
										continue;

									// get the volume and correctly orient the tets
									// TODO: volume for degenerate tets?

									// TODO: why didn't this function work????
									//psi_real vol1 = psi_orient_tet(gpos1, gvel1);
									psi_real adx = gpos1[0].x - gpos1[3].x;
									psi_real bdx = gpos1[1].x - gpos1[3].x;
									psi_real cdx = gpos1[2].x - gpos1[3].x;
									psi_real ady = gpos1[0].y - gpos1[3].y;
									psi_real bdy = gpos1[1].y - gpos1[3].y;
									psi_real cdy = gpos1[2].y - gpos1[3].y;
									psi_real adz = gpos1[0].z - gpos1[3].z;
									psi_real bdz = gpos1[1].z - gpos1[3].z;
									psi_real cdz = gpos1[2].z - gpos1[3].z;
									psi_real vol1 = -0.16666666666666666*(adx*(bdy*cdz-bdz*cdy)
											+bdx*(cdy*adz-cdz*ady)+cdx*(ady*bdz-adz*bdy));
									if(vol1 < 0.0) { // swap two vertices if the vol1ume is negative
										psi_rvec swp = gpos1[0]; gpos1[0] = gpos1[1]; gpos1[1] = swp;
										swp = gvel1[0]; gvel1[0] = gvel1[1]; gvel1[1] = swp;
										vol1 *= -1;
									}
									if(vol1 <= 0.0) continue; // TODO: treat degenerate tets better!!!
									psi_real rho1 = mass1/vol1; 



									
									// clip the tet
									psi_tet_faces_from_verts(tetfaces, gpos1); // TODO: make all of these planes up front
									curpoly = tpoly;
									psi_clip(&curpoly, tetfaces, 4);
									psi_reduce(&curpoly, moments, 1, 0);
									if(moments[0] <= 0.0) continue;

									// set up the vertex weights,
									// depending on the specified weight scheme
									for(v = 0; v < 4; ++v) {
										// mass-weighted for debugging
										// TODO: different factor for the same element??
										vertweights[v] = 0.5*mass0*rho1;
										// inverse r^2, gpos is already relative to the observer position
										//vertweights[v] = tetbuf0.mass[t]/(gpos[v].x*gpos[v].x+gpos[v].y*gpos[v].y+gpos[v].z*gpos[v].z); 
									}
			
									// subtract. Now moments[0->3] contain barycentric moments,
									// properly weighted wrt the original tet vertices
									// add the mass to the grid, linearly interpolating the vertex weights 
									moments[0] -= (moments[1]+moments[2]+moments[3]); 
									for(v = 0; v < 4; ++v) {
										grid->fields[0][p] += vertweights[v]*moments[v]; 
										mtot += vertweights[v]*moments[v]; 
									}		
								}
							}
						}
					}
				}
				rold = rnew;
			}
		}

		if(p%16==0)
			psi_printf("\rPixel %d of %d, %.1f%%", p, grid->n.k, (100.0*p)/grid->n.k);
	}

	psi_rtree_destroy(&rtree);

	printf("Total mass mapped = %.5e\n", mtot);
	printf("Omtot/4pi-1 = %.5e\n", omtot/(FOUR_PI)-1.0);
	printf("Num. failed queries = %d\n", nfail);

	// make sure all elements were used...
	//for(e = 0; e < mesh->nelem; ++e) {
		//if(elemmarks[e]) continue;
		//printf(" Element %d was never used !!!---!!!\n", e);
	//}
	psi_free(elemmarks);

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
