#include "beamtrace.h"

psi_int psi_load_phi_info(psi_metric* metric) {

	psi_int ngrid, nread, ax, i, bpos;
	FILE* f;

	// don't touch it if already loaded
	if(metric->snapnum > 0) 
		return 1;

	// hard-coded from Radek's files for now
	metric->snapnum = 189; 
	metric->n.i = 128;
	metric->n.j = 128;
	metric->n.k = 128;
	for(ax = 0, ngrid = 1; ax < 3; ++ax)
		ngrid *= metric->n.ijk[ax];
	metric->box.x = 0.4;
	metric->box.y = 0.4;
	metric->box.z = 0.4;


    f = fopen("/home/devon/HDD/PS-128-vm/geo/0189", "rb");
	if(!f) return 0;

	metric->phi = malloc(ngrid*sizeof(psi_real));
	metric->gradphi = malloc(3*ngrid*sizeof(psi_real));

	// Radek's data files are float32
	psi_int bufsz = 4096;
	float* buf = malloc(bufsz*sizeof(float));

	// do phi first
	fseek(f, 0, SEEK_SET);
	for(bpos = 0, nread = bufsz; bpos < ngrid; bpos += nread) {
		if(bpos + nread > ngrid)
			nread = ngrid - bpos;
		fread(buf, nread, sizeof(float), f); // read 10 bytes to our buffer
		for(i = 0; i < nread; ++i) 
			metric->phi[bpos+i] = (psi_real)buf[i];
	}


	// grad phi 
	fseek(f, 4*ngrid*sizeof(float), SEEK_SET);
	for(bpos = 0, nread = bufsz; bpos < 3*ngrid; bpos += nread) {
		if(bpos + nread > 3*ngrid)
			nread = 3*ngrid - bpos;
		fread(buf, nread, sizeof(float), f); // read 10 bytes to our buffer
		for(i = 0; i < nread; ++i) 
			metric->gradphi[bpos+i] = (psi_real)buf[i];
	}

	printf("Loaded the metric files!!!\n");

	free(buf);
	fclose(f);
	return 1;
}

void psi_get_phi_and_grad(psi_real* phi, psi_4vec *gradphi, psi_4vec pos, psi_metric* metric) {

	psi_int ax;
	*phi = 0.0;
	memset(gradphi, 0, sizeof(psi_4vec));

	if(metric->type == PSI_METRIC_FLRW) {

		//return;
		if(!psi_load_phi_info(metric)) {
			printf("Failed to load metric!!!\n");
			exit(0);

		}

		psi_rvec xtmp;
		psi_dvec ind;
		for(ax = 0; ax < 3; ++ax) {
			xtmp.xyz[ax] = fmod(pos.txyz[ax+1], metric->box.xyz[ax]);
			if(xtmp.xyz[ax] < 0.0) 
				xtmp.xyz[ax] += metric->box.xyz[ax];
			ind.ijk[ax] = floor(xtmp.xyz[ax]*metric->n.ijk[ax]/metric->box.xyz[ax]);
		}
		psi_int flatind = metric->n.j*metric->n.k*ind.i
			+ metric->n.k*ind.j + ind.k;

		*phi = metric->phi[flatind];
		for(ax = 1; ax < 4; ++ax) 
			gradphi->txyz[ax] = metric->gradphi[3*flatind+(ax-1)]; 

		//*phi *= 1e0; 
		//for(ax = 1; ax < 4; ++ax) 
			//gradphi->txyz[ax] *= 1e0; 

		//printf("Got phi = %.5e, grad = %.5e %.5e %.5e %.5e\n", *phi, gradphi->t, gradphi->x, gradphi->y, gradphi->z);
	
		//psi_real r0 = 2.0;
		//psi_4vec x0 = {{0.0,30.0,20.0,20.0}};
		//psi_real phitmp = -0.1*exp(-((pos.x-x0.x)*(pos.x-x0.x)+(pos.y-x0.y)*(pos.y-x0.y)+(pos.z-x0.z)*(pos.z-x0.z))/(2*r0*r0));
		//*phi = phitmp;
		//for(ax = 1; ax < 4; ++ax) 
			//gradphi->txyz[ax] = phitmp*(x0.txyz[ax]-pos.txyz[ax])/(r0*r0);
	}

}

void psi_get_4acc(psi_4vec *acc, psi_4vec pos, psi_4vec vel, psi_metric* metric) {

	psi_int ax;
	psi_real phi;
	psi_4vec gradphi;

	memset(acc, 0, sizeof(psi_4vec));

	if(metric->type == PSI_METRIC_FLRW) {

		// set Hubble term first
		// TODO: pass this as a metric parameter
		psi_real hubble = 0.0;
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

	else if(metric->type == PSI_METRIC_KERR) {



		psi_real sigma, delta, dsigmadtheta, dsigmadr, ddeltadr;
		psi_real sintheta, costheta, cottheta;
		psi_real sin2theta, cos2theta;
		psi_real M, r, a;
		psi_real ut, ur, utheta, uphi; 

		M = 1;
		a = 0.7;

		// set up preliminary stuff
		ut = vel.t;
		ur = vel.x;
		utheta = vel.y;
		uphi = vel.x;
		r = pos.x;
		sintheta = sin(pos.y);
		costheta = cos(pos.y);
		cottheta = costheta/sintheta;
		printf("cottheta = %.5e\n", cottheta);
		sin2theta = sintheta*sintheta;
	   	cos2theta = costheta*costheta;	
		sigma = r*r + a*a*costheta*costheta;
		delta = r*r - 2*M*r + a*a; 
		dsigmadr = 2*r;
		dsigmadtheta = -2*a*a*costheta*sintheta; 
		ddeltadr = -2*M + 2*r; 

		// this is all Mathematica vomit
	
		// t
		acc->t = (-2*M*(r*utheta*(-2*a*sintheta*sintheta*uphi*sigma*(2*a*a*M*r*sin2theta + (a*a + r*r)*dsigmadtheta) + 
	           ut*(-6*a*a*M*r*sintheta*sintheta*dsigmadtheta + sigma*(8*a*a*M*r*sin2theta + (a*a + r*r)*dsigmadtheta))) + 
	        ur*(ut*(-6*a*a*M*r*sintheta*sintheta + (a*a + r*r)*sigma)*(-sigma + r*dsigmadr) + 
	           2*a*sintheta*sintheta*uphi*sigma*((a*a - r*r)*sigma - r*(a*a + r*r)*dsigmadr))))/
	    (sigma*(12*a*a*M*M*r*r*sintheta*sintheta - M*r*(a*a + 2*r*r + a*a*cos2theta)*sigma + (a*a + r*r)*sigma*sigma));
	
		// r
	   acc->x = (2*M*ut*ut*delta*delta*(sigma - r*dsigmadr) - 8*a*M*sintheta*sintheta*ut*uphi*delta*delta*(sigma - r*dsigmadr) + 
	      2*sintheta*sintheta*uphi*uphi*delta*delta*(a*a*M*sintheta*sintheta*sigma + r*sigma*sigma - a*a*M*r*sintheta*sintheta*dsigmadr) + 
	      sigma*sigma*(2*ur*utheta*(-delta*dsigmadtheta) + utheta*utheta*delta*delta*dsigmadr + 
	         ur*ur*(sigma*ddeltadr - delta*dsigmadr)))/(2.*delta*sigma*sigma*sigma);
	
	   // theta
	   acc->y = (-2*M*r*ut*ut*delta*delta*dsigmadtheta + 8*a*M*r*sintheta*ut*uphi*delta*delta*(-2*costheta*sigma + sintheta*dsigmadtheta) + 
	      uphi*uphi*delta*delta*(8*a*a*M*r*costheta*sintheta*sintheta*sintheta*sigma + (a*a + r*r)*sin2theta*sigma*sigma - 
	         2*a*a*M*r*sintheta*sintheta*sintheta*sintheta*dsigmadtheta) - sigma*sigma*
	       (utheta*utheta*delta*delta*dsigmadtheta + ur*ur*(-delta*dsigmadtheta) + 
	         2*ur*utheta*delta*delta*dsigmadr))/(2.*delta*delta*sigma*sigma*sigma);
	
	   // phi
	   acc->z = (-4*a*M*ut*sigma*(r*utheta*(4*M*r*cottheta - 2*cottheta*sigma + dsigmadtheta) + ur*(-sigma + r*dsigmadr)) + 
	      uphi*(2*utheta*(2*M*r*(r*r + a*a*cos2theta)*cottheta*sigma*sigma - (a*a + r*r)*cottheta*sigma*sigma*sigma + 
	            6*a*a*M*M*r*r*sintheta*sintheta*dsigmadtheta - a*a*M*r*sintheta*sigma*(8*M*r*costheta - sintheta*dsigmadtheta)) + 
	         ur*(M*(-a*a + 4*r*r + a*a*cos2theta)*sigma*sigma - 2*r*sigma*sigma*sigma - 2*a*a*M*r*sintheta*sintheta*sigma*(6*M - dsigmadr) + 
	            12*a*a*M*M*r*r*sintheta*sintheta*dsigmadr)))/
	    (sigma*(12*a*a*M*M*r*r*sintheta*sintheta - M*r*(a*a + 2*r*r + a*a*cos2theta)*sigma + (a*a + r*r)*sigma*sigma));


	}

}


psi_real psi_4dot(psi_4vec dx0, psi_4vec dx1, psi_4vec x, psi_metric* metric) {

	psi_real phi;
	psi_4vec gradphi;

	if(metric->type == PSI_METRIC_MINKOWSKI) {
		return (-dx0.t*dx1.t*PSI_CLIGHT*PSI_CLIGHT + dx0.x*dx1.x + dx0.y*dx1.y + dx0.z*dx1.z);
	}

	if(metric->type == PSI_METRIC_FLRW) {
		psi_real scalefac = 1.0;
		psi_get_phi_and_grad(&phi, &gradphi, x, metric);
		return scalefac*scalefac*(-dx0.t*dx1.t*PSI_CLIGHT*PSI_CLIGHT*(1+2*phi)
				+ (dx0.x*dx1.x + dx0.y*dx1.y + dx0.z*dx1.z)*(1-2*phi));
	}

	if(metric->type == PSI_METRIC_KERR) {


		psi_real sigma, delta;
		psi_real sintheta, costheta;
		psi_real M, r, a;

		M = 1;
		a = 0.7;

		// set up preliminary stuff
		r = x.x;
		sintheta = sin(x.y);
		costheta = cos(x.y);
		sigma = r*r + a*a*costheta*costheta;
		delta = r*r - 2*M*r + a*a; 

		return (dx0.x*dx1.x*sigma*sigma + delta*(2*M*r*(dx0.t*dx1.t - 2*a*(dx0.t*dx1.z + dx0.z*dx1.t)*sintheta*sintheta + a*a*dx0.z*dx1.z*sintheta*sintheta*sintheta*sintheta) + 
				(-(dx0.t*dx1.t) + dx0.z*dx1.z*(a*a + r*r)*sintheta*sintheta)*sigma + dx0.y*dx1.y*sigma*sigma))/(delta*sigma);
	
	}
}

void psi_setup_ray_coords(psi_beam* beam, psi_rvec pos, psi_rvec vel, psi_metric* metric) {

	// TODO: more clever transformation from observer frame!

	psi_int ax;
	psi_real phi;
	psi_4vec gradphi;

	// first, convert Cartesian position coordinates
	// to the appropriate spatial coordinates for this metric
	beam->alpha = 0.0;
	if(metric->type == PSI_METRIC_MINKOWSKI || metric->type == PSI_METRIC_FLRW) {
		// for Cartesian metrics, just use the given spatial coordinates
		for(ax = 0; ax < 3; ++ax) {
			beam->pos.txyz[ax+1] = pos.xyz[ax];
			beam->vel.txyz[ax+1] = vel.xyz[ax];
		}
	}
	else if(metric->type == PSI_METRIC_KERR) {
		// Kerr metric is in Boyer-Lindquist (spherical) coordinates
		// r, theta, phi
		beam->pos.x = sqrt(pos.x*pos.x+pos.y*pos.y+pos.z*pos.z);
		beam->pos.y = acos(pos.z/beam->pos.x); 
		beam->pos.z = atan2(pos.y, pos.x); 

		// TODO: beam 4-velocity from vel!

		// dr/dt
		psi_real x2y2 = pos.x*pos.x + pos.y*pos.y;
		beam->vel.x = (vel.x*pos.x+vel.y*pos.y+vel.z*pos.z)/beam->pos.x;
		beam->vel.x = (vel.x*pos.y-vel.y*pos.x)/x2y2; 
		beam->vel.z = (pos.z*(vel.x*pos.x+vel.y*pos.y)-x2y2*vel.z)/(beam->pos.x*beam->pos.x*sqrt(x2y2));




		////




	}



	// last, set the time coordinate of the 4-velocity
	// to enforce a null trajectory
	beam->pos.t = 0.0; // set observer reference time to zero
	if(metric->type == PSI_METRIC_MINKOWSKI) {
		beam->vel.t = sqrt(beam->vel.x*beam->vel.x 
				+ beam->vel.y*beam->vel.y + beam->vel.z*beam->vel.z)/PSI_CLIGHT;
	}
	else if(metric->type == PSI_METRIC_FLRW) {
		psi_get_phi_and_grad(&phi, &gradphi, beam->pos, metric);
		beam->vel.t = sqrt((beam->vel.x*beam->vel.x + beam->vel.y*beam->vel.y 
				+ beam->vel.z*beam->vel.z)*(1-2*phi)/(1+2*phi))/PSI_CLIGHT;
	}
	else if(metric->type == PSI_METRIC_KERR) {

		psi_real sigma, delta;
		psi_real sintheta, costheta;
		psi_real M, r, a;
		M = 1;
		a = 0.7;

		// set up preliminary stuff
		r = beam->pos.x;
		sintheta = sin(beam->pos.y);
		costheta = cos(beam->pos.y);
		sigma = r*r + a*a*costheta*costheta;
		delta = r*r - 2*M*r + a*a; 
		beam->vel.t = (8*a*beam->vel.z*M*r*sintheta*sintheta*delta + sqrt(64*a*a*beam->vel.z*beam->vel.z*M*M*r*r*sintheta*sintheta*sintheta*sintheta*delta*delta - 
			4*(2*M*r*delta - delta*sigma)*(2*a*a*beam->vel.z*beam->vel.z*M*r*sintheta*sintheta*sintheta*sintheta*delta + a*a*beam->vel.z*beam->vel.z*sintheta*sintheta*delta*sigma + 
			beam->vel.z*beam->vel.z*r*r*sintheta*sintheta*delta*sigma + beam->vel.x*beam->vel.x*sigma*sigma + beam->vel.y*beam->vel.y*delta*sigma*sigma)))/(2.*(2*M*r*delta - delta*sigma));
	}

}

// the main function
// traces beams as defined by the grb_params passed in
void psi_beamtrace(psi_grid* grid, psi_mesh* mesh, psi_int bstep, psi_metric* metric, psi_real* outar, psi_dvec ardim) {

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

		// make the position and velocity 4-vectors
		// enforce ds^2 = 0 for the tangent vector 
		psi_setup_ray_coords(&cbeam, obspos, center, metric);

		// trace the beam
		dalpha = 0.01;
		for(terminate_beam = 0, step = 0; step < 100000 && !terminate_beam; ++step) {

			// save ray samples for plotting
			if(metric->type == PSI_METRIC_KERR) {
				psi_real r = cbeam.pos.x;
				psi_real sintheta = sin(cbeam.pos.y);
				psi_real costheta = cos(cbeam.pos.y);
				psi_real sinphi = sin(cbeam.pos.z);
				psi_real cosphi = cos(cbeam.pos.z);
				outar[p*ardim.j*ardim.k + step*ardim.k + 0] = cbeam.pos.t; // t
				outar[p*ardim.j*ardim.k + step*ardim.k + 1] = r*sintheta*cosphi; 
				outar[p*ardim.j*ardim.k + step*ardim.k + 2] = r*sintheta*sinphi;
				outar[p*ardim.j*ardim.k + step*ardim.k + 3] = r*costheta;

			}
			else {
				for(ax = 0; ax < 4; ++ax)
					outar[p*ardim.j*ardim.k + step*ardim.k + ax] = cbeam.pos.txyz[ax];
			}

			if(isnan(outar[p*ardim.j*ardim.k + step*ardim.k + 0])) {
				break;
			}

			//printf("step %d, Metric = %d, pos = %f %f %f %f, vel = %f %f %f %f\n", step, metric->type, outar[p*ardim.j*ardim.k + step*ardim.k + 0],outar[p*ardim.j*ardim.k + step*ardim.k + 1],outar[p*ardim.j*ardim.k + step*ardim.k + 2],outar[p*ardim.j*ardim.k + step*ardim.k + 3], cbeam.vel.t,cbeam.vel.x,cbeam.vel.y,cbeam.vel.z);
			//printf(" - alpha = %.5e, ds2 = %.5e\n", cbeam.alpha, psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric));
			////if(step > 300)
				////exit(0);
		
		
			// get the 4-acceleration
			// calculate the displacement vector for this step
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
			if(cbeam.alpha > 100.0)
				terminate_beam = 1;
		}

		// check the final 4-velocity to see that it is still null
		//printf("Final velocity U*u = %.5e\n", psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric));
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

