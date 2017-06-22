#include "skymap.h"


// the main function
// traces beams as defined by the grb_params passed in
void psi_skymap(psi_grid* grid, psi_int bstep) {

	psi_int p, b, ax; 
	psi_dvec grind;
	psi_real vol;
	psi_rvec corners[64], center;
	psi_beam beam;

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
				beam.pos[0][ax+1] = corners[b].xyz[ax];
				beam.pos[1][ax+1] = corners[(b+1)%(4*bstep)].xyz[ax];
				beam.pos[2][ax+1] = center.xyz[ax]; 
			}
	
			// TODO: rotate/transform/boost the beam position to the source
	
			// trace the beam, then accumulate the
			// resulting flux into the pixel
			//trace_beam(grbp, &beam);
		}
		beam.flux = vol;
		grid->fields[0][p] += beam.flux; 
	}

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
