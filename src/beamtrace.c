#include "beamtrace.h"

void psi_get_phi_and_grad(psi_real* phi, psi_4vec *gradphi, psi_4vec pos, psi_int metric) {

	psi_int ax;
	*phi = 0.0;
	memset(gradphi, 0, sizeof(psi_4vec));

	if(metric == PSI_METRIC_FLRW) {

		//return;
	
		psi_real r0 = 2.0;
		psi_4vec x0 = {{0.0,30.0,20.0,20.0}};
		psi_real phitmp = -0.2*exp(-((pos.x-x0.x)*(pos.x-x0.x)+(pos.y-x0.y)*(pos.y-x0.y)+(pos.z-x0.z)*(pos.z-x0.z))/(2*r0*r0));
		*phi = phitmp;
		for(ax = 1; ax < 4; ++ax) 
			gradphi->txyz[ax] = phitmp*(x0.txyz[ax]-pos.txyz[ax])/(r0*r0);
	}

}

void psi_get_4acc(psi_4vec *acc, psi_4vec pos, psi_4vec vel, psi_int metric) {

	psi_int ax;
	psi_real phi;
	psi_4vec gradphi;

	if(metric == PSI_METRIC_MINKOWSKI) {
		for(ax = 0; ax < 4; ++ax)
			acc->txyz[ax] = 0.0;
	}


	else if(metric == PSI_METRIC_FLRW) {

		// Hubble parameter
		psi_real hubble = 0.0;

		// set Hubble term first
		// TODO: the t-part of this has been reduce using the null condition
		// -- should we not enforce it here??
		for(ax = 0; ax < 4; ++ax)
			acc->txyz[ax] = -2*hubble*vel.txyz[ax]*vel.t; 

		// do the gradients of phi
		psi_get_phi_and_grad(&phi, &gradphi, pos, metric);
		psi_real vel2 = vel.t*vel.t + vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
		psi_real vdotgrad = vel.t*gradphi.t + vel.x*gradphi.x + vel.y*gradphi.y + vel.z*gradphi.z;
		acc->t += -1/(1+2*phi)*(2*vel.t*vdotgrad-vel2*gradphi.t);
		for(ax = 1; ax < 4; ++ax)
			acc->txyz[ax] += 1/(1-2*phi)*(2*vel.txyz[ax]*vdotgrad-vel2*gradphi.txyz[ax]);

	}

}


psi_real psi_4dot(psi_4vec dx0, psi_4vec dx1, psi_4vec x, psi_int metric) {

	psi_real phi;
	psi_4vec gradphi;

	if(metric == PSI_METRIC_MINKOWSKI) {
		return (-dx0.t*dx1.t*PSI_CLIGHT*PSI_CLIGHT + dx0.x*dx1.x + dx0.y*dx1.y + dx0.z*dx1.z);
	}

	if(metric == PSI_METRIC_FLRW) {
		psi_real scalefac = 1.0;
		psi_get_phi_and_grad(&phi, &gradphi, x, metric);
		return scalefac*scalefac*(-dx0.t*dx1.t*PSI_CLIGHT*PSI_CLIGHT*(1+2*phi)
				+ (dx0.x*dx1.x + dx0.y*dx1.y + dx0.z*dx1.z)*(1-2*phi));
	}
}

void psi_enforce_null_time(psi_4vec *vel, psi_4vec x, psi_int metric) {

	// TODO: more clever transformation from observer frame!

	psi_real phi;
	psi_4vec gradphi;

	if(metric == PSI_METRIC_MINKOWSKI) {
		vel->t = sqrt(vel->x*vel->x + vel->y*vel->y + vel->z*vel->z)/PSI_CLIGHT;
	}

	else if(metric == PSI_METRIC_FLRW) {
		psi_get_phi_and_grad(&phi, &gradphi, x, metric);
		vel->t = sqrt((vel->x*vel->x + vel->y*vel->y + vel->z*vel->z)*(1-2*phi)/(1+2*phi))/PSI_CLIGHT;
	}


}

// the main function
// traces beams as defined by the grb_params passed in
void psi_beamtrace(psi_grid* grid, psi_mesh* mesh, psi_int bstep, psi_int metric, psi_real* outar, psi_dvec ardim) {

	psi_int p, b, ax, e, v, tind, t, terminate_beam, step; 
	psi_dvec grind;
	psi_real dalpha, ds2tot, vol;
	psi_real moments[10];
	psi_real vertweights[4];
	psi_rvec corners[64], center, brbox[2];
	psi_rvec gpos[mesh->dim+1], grbox[2];
	psi_int vpere = mesh->elemtype;
	psi_beam beam;
	psi_plane beamfaces[16];
	psi_tet_buffer tetbuf;
	psi_beam cbeam;

	psi_rvec rawverts[3]; 
	psi_rvec beamverts[6];

	psi_rvec tpos[vpere], tvel[vpere], trbox[2];
	psi_real tmass;

	psi_4vec dx, dv, acc;

	if(grid->type != PSI_GRID_HPRING)
		return;

	psi_rvec obspos = {{20.,20.,20.}};


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

		if(fabs(center.z) > 0.1) continue;
		//if(p > 10) continue;
		//printf("pix %d, center = %f %f %f\n", p, center.x, center.y, center.z);


		// make the position and velocity 4-vectors
		// enforce ds^2 = 0 for the tangent vector 
		for(ax = 0; ax < 3; ++ax) {
			cbeam.pos.txyz[ax+1] = obspos.xyz[ax];
			cbeam.vel.txyz[ax+1] = center.xyz[ax];
		}
		cbeam.pos.t = 0.0;
		psi_enforce_null_time(&cbeam.vel, cbeam.pos, metric);

		// trace the beam
		dalpha = 0.01;
		//psi_real ds2tot = 0.0;
		cbeam.alpha = 0.0;
		for(terminate_beam = 0, step = 0; step < 100000 && !terminate_beam; ++step) {
		
			// save ray samples for plotting
			for(ax = 0; ax < 4; ++ax)
				outar[p*ardim.j*ardim.k + step*ardim.k + ax] = cbeam.pos.txyz[ax];
		
			// get the 4-acceleration
			// calculate the displacement vector for this step
			// get ds2 to check enforcement of null condition
			psi_get_4acc(&acc, cbeam.pos, cbeam.vel, metric);
			for(ax = 0; ax < 4; ++ax) {
				dx.txyz[ax] = dalpha*cbeam.vel.txyz[ax]; 
				dv.txyz[ax] = dalpha*acc.txyz[ax]; 
			} 

			// update beam 4-vectors 
			// TODO: forward Euler for now
			cbeam.alpha += dalpha;
			for(ax = 0; ax < 4; ++ax) {
				cbeam.pos.txyz[ax] += dx.txyz[ax];
				cbeam.vel.txyz[ax] += dv.txyz[ax];
			} 
		
			// check beam stopping criteria
			if(cbeam.alpha > 30.0)
				terminate_beam = 1;
		}

		// check the final 4-velocity to see that it is still null
		printf("Final velocity U*u = %.5e\n", psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric));
		grid->fields[0][p] = psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric);

#if 0

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
			psi_real domega = psi_omega3(rawverts[0], rawverts[1], rawverts[2]);
			if(domega <= 0.0)
				printf("Bad domega!!!!\n");

			omtot += domega;
	
			psi_real rmin = 15.0;
			psi_real rmax = 15.0 + 0.001;
			psi_real rstep = .001001; 
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

					psi_real ttmass = 0.0;
					elemmarks[e] = 1;

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

						ttmass += tetbuf.mass[t]; 
			
						// copy the position to the ghost array and compute its aabb
						memcpy(gpos, &tetbuf.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
						psi_aabb(gpos, mesh->dim+1, grbox);

						if(mode == PSI_SKYMAP_RHO_LINEAR) {

							// TODO: check the boxes before clipping....
							psi_clip_reduce_tet(gpos, beamfaces, 5, moments);
							if(moments[0] <= 0.0) {
								//if(moments[0] < 0.0)
									//psi_printf("Moments < 0 ! %.5e\n", moments[0]);
								continue;
							} 
	
	
							// set up the vertex weights,
							// depending on the specified weight scheme
							for(v = 0; v < 4; ++v) {
	
								// mass-weighted for debugging
								vertweights[v] = tetbuf.mass[t]; 
	
								// inverse r^2, gpos is already relative to the observer position
								//vertweights[v] = tetbuf.mass[t]/(gpos[v].x*gpos[v].x+gpos[v].y*gpos[v].y+gpos[v].z*gpos[v].z); 
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


						else if(mode == PSI_SKYMAP_RHO_SQUARED) {

							//printf(" tet grbox = %f %f %f to %f %f %f \n", grbox[0].x, grbox[0].y, grbox[0].z, grbox[1].x, grbox[1].y, grbox[1].z);

							// make a new query to find elements close to this one
							psi_int e1;
							psi_rtree_query qry1;
							psi_rtree_query_init(&qry1, &rtree, grbox);
							while(psi_rtree_query_next(&qry1, &e1)) {

								printf("Found possibly overlapping elements!\n");


							}


#if 0

							// TODO: check the boxes before clipping....
							psi_clip_reduce_tet(gpos, beamfaces, 5, moments);
							if(moments[0] <= 0.0) {
								//if(moments[0] < 0.0)
									//psi_printf("Moments < 0 ! %.5e\n", moments[0]);
								continue;
							} 
	
	
							// set up the vertex weights,
							// depending on the specified weight scheme
							for(v = 0; v < 4; ++v) {
	
								// mass-weighted for debugging
								vertweights[v] = tetbuf.mass[t]; 
	
								// inverse r^2, gpos is already relative to the observer position
								//vertweights[v] = tetbuf.mass[t]/(gpos[v].x*gpos[v].x+gpos[v].y*gpos[v].y+gpos[v].z*gpos[v].z); 
							}
	
							// subtract. Now moments[0->3] contain barycentric moments,
							// properly weighted wrt the original tet vertices
							// add the mass to the grid, linearly interpolating the vertex weights 
							moments[0] -= (moments[1]+moments[2]+moments[3]); 
							for(v = 0; v < 4; ++v) {
								grid->fields[0][p] += vertweights[v]*moments[v]; 
								mtot += vertweights[v]*moments[v]; 
							}		
#endif			
						
						}



					}
				}
				rold = rnew;
			}
		}
#endif

		if(p%512==0)
			psi_printf("\rPixel %d of %d, %.1f%%", p, grid->n.k, (100.0*p)/grid->n.k);

#ifdef PYMODULE 
		PyErr_CheckSignals();
#endif
	}

	// if we are in healpix, solve for the correction factors
	// that make it truly equal-area
	//if(grbp->pix_type == PIX_HEALPIX)
		//correct_healpix_edges(grbp, 1.0e-6);
} 

