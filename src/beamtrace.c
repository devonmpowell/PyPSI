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

void psi_get_updated_4vel(psi_4vec *newvel, psi_4vec pos, psi_4vec vel, psi_real dalpha, psi_metric* metric) {

	psi_int ax;
	psi_real phi;
	psi_4vec gradphi, acc;

	memset(newvel, 0, sizeof(psi_4vec));

	if(metric->type == PSI_METRIC_FLRW) {

		// set Hubble term first
		// TODO: pass this as a metric parameter
		psi_real hubble = 0.0;
		for(ax = 0; ax < 4; ++ax)
			acc.txyz[ax] = -2*hubble*vel.txyz[ax]*vel.t; 

		// do the gradients of phi
		psi_get_phi_and_grad(&phi, &gradphi, pos, metric);
		psi_real vel2 = vel.t*vel.t + vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
		psi_real vdotgrad = vel.t*gradphi.t + vel.x*gradphi.x + vel.y*gradphi.y + vel.z*gradphi.z;
		acc.t += -1/(1+2*phi)*(2*vel.t*vdotgrad-vel2*gradphi.t);
		for(ax = 1; ax < 4; ++ax)
			acc.txyz[ax] += 1/(1-2*phi)*(2*vel.txyz[ax]*vdotgrad-vel2*gradphi.txyz[ax]);

		// update the 4-velocity
		*newvel = vel;
		for(ax = 0; ax < 4; ++ax)
			newvel->txyz[ax] += dalpha*acc.txyz[ax]; 

	}

	else if(metric->type == PSI_METRIC_KERR) {

		psi_real M = 1.0;
		psi_real a = 0.7;

		// polar coordinates
		psi_real r = pos.x; 
		psi_real theta = pos.y; 
		psi_real rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		psi_real delta = r*r - 2*r + a*a;

		// constants of motion
		psi_real k = vel.x;
		psi_real h = vel.y;
		psi_real Q = vel.z;



		// tdot
		psi_real pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

		// phidot
		psi_real pphi = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
		pphi = pphi / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

		// thetadot
		psi_real thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
		thetadotsq = thetadotsq / (rhosq*rhosq);


		psi_int thetadot_sign = 1;
		psi_int rdot_sign = 1;
		//if(thetadotsq < 0 && thetasign_count >= COUNT_MIN)
		//{
			//thetadot_sign *= -1;
			//thetasign_count = 0;
			//return;
		//}
		//if (thetasign_count <= COUNT_MIN) thetasign_count++;

		// take the square roots and get the right signs
		psi_real ptheta = sqrtf(abs(thetadotsq)) * thetadot_sign;

		// rdot
		psi_real rdotsq = k*pt - h*pphi - rhosq*ptheta*ptheta;
		rdotsq = rdotsq * delta/rhosq;

		//if(rdotsq < 0 && rsign_count >= COUNT_MIN)
		//{
			//rdot_sign *= -1;
			//rsign_count = 0;
			//return;
		//}
		//if (rsign_count <= COUNT_MIN) rsign_count++;

		psi_real pr = sqrtf(abs(rdotsq)) * rdot_sign;

#if 0
		step = abs( (r-(T)horizon)/pr ) / tol;
		// if the step is smaller in theta or phi (near 0/pi/2pi), use that instead
		if( step > abs( theta/ptheta ) / tol ) step = abs( theta/ptheta ) / tol;
//		if( step > abs( phi/pphi ) / tol ) step = abs( phi/pphi ) / tol;
//		if( step > abs( (phi - M_PI)/pphi ) / tol ) step = abs( (phi - M_PI)/pphi ) / tol;
//		if( step > abs( (phi - 2*M_PI)/pphi ) / tol ) step = abs( (phi - 2*M_PI)/pphi ) / tol;
		// don't let the step be stupidly small
		if( step < MIN_STEP ) step = MIN_STEP;
		if(reverse) step *= -1;
#endif

		// calculate new position
#if 0
		if(r <= horizon) break;
		if(boundary > 0 && r <= boundary) break;
#endif

		newvel->t = pt;
		newvel->x = pr;
		newvel->y = ptheta;
		newvel->z = pphi;


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

		// TODO: correct units!
		return (dx0.x*dx1.x*sigma*sigma + delta*(2*M*r*(dx0.t*dx1.t - 2*a*(dx0.t*dx1.z + dx0.z*dx1.t)*sintheta*sintheta + a*a*dx0.z*dx1.z*sintheta*sintheta*sintheta*sintheta) + 
				(-(dx0.t*dx1.t) + dx0.z*dx1.z*(a*a + r*r)*sintheta*sintheta)*sigma + dx0.y*dx1.y*sigma*sigma))/(delta*sigma);
	
	}
}

void psi_setup_ray_coords(psi_4vec* pos, psi_4vec* vel, psi_int* info, psi_rvec cartpos, psi_rvec cartvel, psi_metric* metric) {

	// TODO: more clever transformation from observer frame!

	psi_int ax;
	psi_real phi;
	psi_4vec gradphi;

	// first, convert Cartesian position coordinates
	// to the appropriate spatial coordinates for this metric
	if(metric->type == PSI_METRIC_MINKOWSKI || metric->type == PSI_METRIC_FLRW) {
		// for Cartesian metrics, just use the given Cartesian coordinates
		for(ax = 0; ax < 3; ++ax) {
			pos->txyz[ax+1] = cartpos.xyz[ax];
			vel->txyz[ax+1] = cartvel.xyz[ax];
		}
	}
	else if(metric->type == PSI_METRIC_KERR) {
		// Kerr metric is in Boyer-Lindquist (spherical) coordinates
		// r, theta, phi

		psi_real M = 1.0;
		psi_real a = 0.7;

		psi_real V = 0.0;

		// the beam spherical coordinates
		psi_real r = sqrt(cartpos.x*cartpos.x+cartpos.y*cartpos.y+cartpos.z*cartpos.z);
		psi_real theta = acos(cartpos.z/cartpos.x); 
		psi_real phi = atan2(cartpos.y, cartpos.x); 

		// Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
		// of a source at (t,r,theta,phi) orbiting azimuthally at angular velocity V
		psi_real rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		psi_real delta = r*r - 2*r + a*a;
		psi_real sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);
	
		// metric coefficients
		psi_real e2nu = rhosq * delta / sigmasq;
		psi_real e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
		psi_real omega = 2*a*r / sigmasq;

		//printf(" e2nu, e2psi, omega = %.5e %.5e %.5e\n", e2nu, e2psi, omega); 
	
		// tetrad basis vector components
		psi_real et0 = (1/sqrtf(e2nu))/sqrtf(1 - (V - omega)*(V - omega)*e2psi/e2nu);
		psi_real et3 = (1/sqrtf(e2nu))*V / sqrtf(1 - (V - omega)*(V - omega)*e2psi/e2nu);
		psi_real e10 = (V - omega)*sqrtf(e2psi/e2nu) / sqrtf(e2nu - (V-omega)*(V-omega)*e2psi);
		psi_real e13 = (1/sqrtf(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrtf(e2nu - (V-omega)*(V-omega)*e2psi);
		psi_real e22 = 1/sqrtf(rhosq);
		psi_real e31 = sqrtf(delta/rhosq);
	
		// photon 4-momentum in source frame
		psi_real E = 1.0;
		//psi_real rdotprime[] = { E, E*sin(alpha)*sin(beta), E*sin(alpha)*cos(beta), E*cos(alpha) };
		psi_real rdotprime[] = { E, E*cartvel.x, E*cartvel.y, E*cartvel.z};
	
		psi_real tdot = rdotprime[0]*et0 + rdotprime[1]*e10;
		psi_real phidot = rdotprime[0]*et3 + rdotprime[1]*e13;
		psi_real rdot = rdotprime[3]*e31;
		psi_real thetadot = rdotprime[2]*e22;


		//printf(" et3, e13 = %f %f\n", et3, e13); 
		//printf(" tdot, phidot, thetadot = %f %f %f\n", tdot, phidot, thetadot);

	
		// find the corresponding values of k, h and Q using the geodesic equations
		psi_real k = (1 - 2*r/rhosq)*tdot + (2*a*r*sin(theta)*sin(theta)/rhosq)*phidot;
		psi_real h = phidot * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
		h = h - 2*a*r*k*sin(theta)*sin(theta);
		h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );
		psi_real Q = rhosq*rhosq*thetadot*thetadot - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));

		// set the beam coordinates in BL system
		pos->x = r;
		pos->y = theta;
		pos->z = phi;

		vel->x = k;
		vel->y = h;
		vel->z = Q;

#define sgn(x) (((x) > 0.0)? 1: -1)
		info[0] = 0;
		info[1] = sgn(rdot);
		info[2] = 0;
		info[3] = sgn(thetadot);

		printf(" pos = %f %f %f , khQ = %f %f %f\n", pos->x, pos->y, pos->z, vel->x, vel->y, vel->z);

	}



	// last, set the time coordinate of the 4-velocity
	// to enforce a null trajectory
	pos->t = 0.0; // set observer reference time to zero
	if(metric->type == PSI_METRIC_MINKOWSKI) {
		vel->t = sqrt(vel->x*vel->x 
				+ vel->y*vel->y + vel->z*vel->z)/PSI_CLIGHT;
	}
	else if(metric->type == PSI_METRIC_FLRW) {
		psi_get_phi_and_grad(&phi, &gradphi, *pos, metric);
		vel->t = sqrt((vel->x*vel->x + vel->y*vel->y 
				+ vel->z*vel->z)*(1-2*phi)/(1+2*phi))/PSI_CLIGHT;
	}
	else if(metric->type == PSI_METRIC_KERR) {

		psi_real sigma, delta;
		psi_real sintheta, costheta;
		psi_real M, r, a;
		M = 1;
		a = 0.7;

		// set up preliminary stuff
		r = pos->x;
		sintheta = sin(pos->y);
		costheta = cos(pos->y);
		sigma = r*r + a*a*costheta*costheta;
		delta = r*r - 2*M*r + a*a; 
		vel->t = (8*a*vel->z*M*r*sintheta*sintheta*delta + sqrt(64*a*a*vel->z*vel->z*M*M*r*r*sintheta*sintheta*sintheta*sintheta*delta*delta - 
			4*(2*M*r*delta - delta*sigma)*(2*a*a*vel->z*vel->z*M*r*sintheta*sintheta*sintheta*sintheta*delta + a*a*vel->z*vel->z*sintheta*sintheta*delta*sigma + 
			vel->z*vel->z*r*r*sintheta*sintheta*delta*sigma + vel->x*vel->x*sigma*sigma + vel->y*vel->y*delta*sigma*sigma)))/(2.*(2*M*r*delta - delta*sigma));
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

	psi_4vec dx, dv, acc, vel;

	if(grid->type != PSI_GRID_HPRING)
		return;

	psi_rvec obspos = {{10.,9.,8.}};


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

		// initialize the beam,
		// make the position and velocity 4-vectors
		// in the correct coordinates for the given metric
		memset(&cbeam, 0, sizeof(psi_beam));
		for(v = 0; v < 3; ++v)
			psi_setup_ray_coords(&cbeam.verts[v].pos, &cbeam.verts[v].vel, cbeam.verts[v].info, obspos, center, metric);

		for(v = 0; v < 3; ++v)
			printf("setup 4-velocity, 4-pos = %f %f %f %f, 4 vel = %f %f %f %f \n", p, cbeam.verts[v].pos.t, cbeam.verts[v].pos.x, cbeam.verts[v].pos.y, cbeam.verts[v].pos.z, cbeam.verts[v].vel.t, cbeam.verts[v].vel.x, cbeam.verts[v].vel.y, cbeam.verts[v].vel.z);

		// trace the beam
		dalpha = 0.01;
		for(terminate_beam = 0, step = 0; step < 1000 && !terminate_beam; ++step) {

			//printf(" pixel %d, step %d\n", p, step);

			// save ray samples for plotting
			psi_4vec cpos;
			for(ax = 0; ax < 4; ++ax) {
				cpos.txyz[ax] = 0.0;
				for(v = 0; v < 3; ++v)
					cpos.txyz[ax] += 0.333333333333*cbeam.verts[v].pos.txyz[ax];
			}
			if(metric->type == PSI_METRIC_KERR) {
				psi_real r = cpos.x;
				psi_real sintheta = sin(cpos.y);
				psi_real costheta = cos(cpos.y);
				psi_real sinphi = sin(cpos.z);
				psi_real cosphi = cos(cpos.z);
				outar[p*ardim.j*ardim.k + step*ardim.k + 0] = cpos.t; // t
				outar[p*ardim.j*ardim.k + step*ardim.k + 1] = r*sintheta*cosphi; 
				outar[p*ardim.j*ardim.k + step*ardim.k + 2] = r*sintheta*sinphi;
				outar[p*ardim.j*ardim.k + step*ardim.k + 3] = r*costheta;

			}
			else {
				for(ax = 0; ax < 4; ++ax)
					outar[p*ardim.j*ardim.k + step*ardim.k + ax] = cpos.txyz[ax];
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

			// for each vertex of the beam...
			for(v = 0; v < 3; ++v) {
			
				psi_get_updated_4vel(&vel, cbeam.verts[v].pos, cbeam.verts[v].vel, dalpha, metric);
				//printf("vert %d 4-vel = %f %f %f %f\n", v, vel.t, vel.x, vel.y, vel.z);

				//printf("pixel %d, 4-pos = %f %f %f %f, 4 vel = %f %f %f %f \n", p, cbeam.verts[v].pos.t, cbeam.verts[v].pos.x, cbeam.verts[v].pos.y, cbeam.verts[v].pos.z, cbeam.verts[v].vel.t, cbeam.verts[v].vel.x, cbeam.verts[v].vel.y, cbeam.verts[v].vel.z);

				// forward Euler 
				//cbeam.alpha += dalpha;
				for(ax = 0; ax < 4; ++ax) {
					dx.txyz[ax] = dalpha*vel.txyz[ax]; //cbeam.vel.txyz[ax]; 
				} 

				// update beam 4-vectors 
				// TODO: forward Euler for now
				for(ax = 0; ax < 4; ++ax) {
					cbeam.verts[v].pos.txyz[ax] += dx.txyz[ax];
					//cbeam.verts[v].vel.txyz[ax] += dv.txyz[ax];
				} 

			
			
			}
		
			// check beam stopping criteria
			//if(cbeam.alpha > 100.0)
				//terminate_beam = 1;
		}

		// check the final 4-velocity to see that it is still null
		//printf("Final velocity U*u = %.5e\n", psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric));
		grid->fields[0][p] = 1.0;//psi_4dot(cbeam.vel, cbeam.vel, cbeam.pos, metric);

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
		if(PyErr_CheckSignals() < 0)
			return;
#endif
	}

	// if we are in healpix, solve for the correction factors
	// that make it truly equal-area
	//if(grbp->pix_type == PIX_HEALPIX)
		//correct_healpix_edges(grbp, 1.0e-6);
} 

