#ifdef HAVE_FFTW
#include "fft.h"

// Hybrid eigensolver routines from J. Kopp (2008)
// http://arxiv.org/abs/physics/0610206
//extern "C" {
	////#include "dsyevh3.h"
	////#include "dsyevc3.h"
	//#include "dsyevj3.h"
//}
#define ind3(i,j,k) (dims.j*dims.k*((i+dims.i)%dims.i) + dims.k*((j+dims.j)%dims.j) + (k+dims.k)%dims.k)

// forward declarations
//void read_hdf5(string filename, string fieldname);
//void read_fortran(string filename);
//void read_gadget_header(string filename);
//void do_phi();
//void do_minmax();

// metadata
psi_real BoxSize; // In Mpc/h
psi_real ScaleFac; // scale factor a
psi_real HubbleParam; // H0 (present-day)
psi_real Omega_m; // Omega_m (present-day)
psi_real Omega_L; // Omega_Lambda (present-day)

#define SMOOTH 0.0

#if 0
// global data arrays
int dims.i, dims.j, dims.k;
psi_real dx, dy, dz;
vector<myreal> rho;

// potential and Hessian eigenvalues 
vector<myreal> phi;
vector<myreal> e0;
vector<myreal> e1;
vector<myreal> e2;

// local minima/maxima
vector<int> minima_phi;
vector<int> maxima_phi;
vector<int> minima_det;
vector<int> maxima_det;

	
	// Correct for hubble expansion.
	// TODO: Is this right?
	psi_real norm = 1.5 * HubbleParam * HubbleParam * Omega_m * (Omega_m/(ScaleFac*ScaleFac*ScaleFac) + Omega_L) /(ScaleFac*ScaleFac*ScaleFac);

	for(int i = 0; i < dims.i*dims.j*dims.k; ++i) {
		rho[i] *= norm/rhobar;
	}
#endif
	
// Gaussian smoothing kernel
// sigma in units of box size
inline psi_real smooth(psi_rvec k, psi_real sigma) {
	//return exp(-0.5*(k.x*k.x + k.y*k.y + k.z*k.z)*sigma*sigma); 
	return 1.0;
}

// Green's fn. for the 13-point Laplacian stencil, to preserve the trace
// of the Hessian to machine precision
inline psi_real laplace(psi_rvec k, psi_rvec d) {
	psi_real fx, fy, fz;
	fx = 1.0/(12.0*d.x*d.x)*(-2*cos(2*k.x*d.x)+32*cos(k.x*d.x)-30);
	fy = 1.0/(12.0*d.y*d.y)*(-2*cos(2*k.y*d.y)+32*cos(k.y*d.y)-30);
	fz = 1.0/(12.0*d.z*d.z)*(-2*cos(2*k.z*d.z)+32*cos(k.z*d.z)-30);
	return 1.0/(fx+fy+fz);
}

void psi_do_phi(psi_grid* grid, psi_real* phi_out, psi_real Gn) {

	fftw_plan p, pinv;
	fftw_complex* rhok;
	psi_int ax;
	psi_dvec halfn, dims;
	psi_rvec kvec, L, dx;
	psi_real kvec2;

	for(ax = 0; ax < 3; ++ax) {
		dims.ijk[ax] = grid->n.ijk[ax];
		halfn.ijk[ax] = dims.ijk[ax]/2 + 1;
		L.xyz[ax] = grid->window[1].xyz[ax]-grid->window[0].xyz[ax];
		dx.xyz[ax] = L.xyz[ax]/dims.ijk[ax];
	}
	
	rhok = (fftw_complex*) fftw_malloc(dims.i*dims.j*halfn.k*sizeof(fftw_complex));
	p = fftw_plan_dft_r2c_3d(dims.i, dims.j, dims.k, grid->fields[0], rhok, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	// Phi
	for(int i = 0; i < dims.i; ++i) 
	for(int j = 0; j < dims.j; ++j) 
	for(int k = 0; k < halfn.k; ++k) {

		// zero the DC component
		if(i == 0 && j == 0 && k == 0) {
			rhok[0][0] = 0.0;
			rhok[0][1] = 0.0;
			continue;
		}
	
		// compute wavenumbers (arbitrary units)	
		kvec.x = TWO_PI/L.x * ((i < halfn.i)? i : (i - dims.i));
		kvec.y = TWO_PI/L.y * ((j < halfn.j)? j : (j - dims.j));
		kvec.z = TWO_PI/L.z * k;
		kvec2 = kvec.x*kvec.x + kvec.y*kvec.y + kvec.z*kvec.z;

		// Filter the FFT
		rhok[dims.j*halfn.k*i + halfn.k*j + k][0] *= FOUR_PI*Gn/(dims.i*dims.j*dims.k)*laplace(kvec, dx)*smooth(kvec, SMOOTH);
		rhok[dims.j*halfn.k*i + halfn.k*j + k][1] *= FOUR_PI*Gn/(dims.i*dims.j*dims.k)*laplace(kvec, dx)*smooth(kvec, SMOOTH);

	}
	pinv = fftw_plan_dft_c2r_3d(dims.i, dims.j, dims.k, rhok, phi_out, FFTW_ESTIMATE);
	fftw_execute(pinv);
	fftw_destroy_plan(pinv);
	fftw_free(rhok);

#if 0
	printf("  Finite difference and eigenvalues...\n");
	
	psi_real phi_xx, phi_yy, phi_zz, phi_xy, phi_xz, phi_yz;
	
	e0.reserve(dims.i*dims.j*dims.k); 
	e1.reserve(dims.i*dims.j*dims.k); 
	e2.reserve(dims.i*dims.j*dims.k); 
	psi_real T[3][3]; 
	psi_real Q[3][3]; // eigenvectors, simply a dummy for now 
	psi_real w[3];

	for(int i = 0; i < dims.i; ++i) 
	for(int j = 0; j < dims.j; ++j) 
	for(int k = 0; k < dims.k; ++k) {
				
		// 2nd-order FD Hessian, from Abramowitz & Stegun
		phi_xx = (-phi[ind3(i+2,j,k)] + 16*phi[ind3(i+1,j,k)] - 30*phi[ind3(i,j,k)] + 16*phi[ind3(i-1,j,k)] - phi[ind3(i-2,j,k)])/(12*dx*dx);
		phi_yy = (-phi[ind3(i,j+2,k)] + 16*phi[ind3(i,j+1,k)] - 30*phi[ind3(i,j,k)] + 16*phi[ind3(i,j-1,k)] - phi[ind3(i,j-2,k)])/(12*dy*dy);
		phi_zz = (-phi[ind3(i,j,k+2)] + 16*phi[ind3(i,j,k+1)] - 30*phi[ind3(i,j,k)] + 16*phi[ind3(i,j,k-1)] - phi[ind3(i,j,k-2)])/(12*dz*dz);
		phi_xy = (phi[ind3(i+1,j+1,k)] - phi[ind3(i+1,j-1,k)] - phi[ind3(i-1,j+1,k)] + phi[ind3(i-1,j-1,k)])/(4*dx*dy);
		phi_xz = (phi[ind3(i+1,j,k+1)] - phi[ind3(i+1,j,k-1)] - phi[ind3(i-1,j,k+1)] + phi[ind3(i-1,j,k-1)])/(4*dx*dz);
		phi_yz = (phi[ind3(i,j+1,k+1)] - phi[ind3(i,j+1,k-1)] - phi[ind3(i,j-1,k+1)] + phi[ind3(i,j-1,k-1)])/(4*dy*dz);

		// This is the COMOVING tidal tensor!
		T[0][0] = phi_xx;
		T[1][1] = phi_yy;
		T[2][2] = phi_zz;
		T[0][1] = phi_xy;
		T[1][0] = phi_xy;
		T[0][2] = phi_xz;
		T[2][0] = phi_xz;
		T[1][2] = phi_yz;
		T[2][1] = phi_yz;

#ifdef PHYSICALHESSIAN
		// rescale comoving Hessian into physical units of length
		for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 3; ++jj) {
			T[ii][jj] /= ScaleFac*ScaleFac*ScaleFac;
		}
		// Modify the trace appropriately to correctly treat Hubble expansion
		// 1/a * (d^2 a)/(dt^2)
		//psi_real addoa = -HubbleParam*HubbleParam*(0.5*Omega_m - Omega_L); 
		psi_real addoa = -HubbleParam*HubbleParam*(0.5*Omega_m/(ScaleFac*ScaleFac*ScaleFac) - Omega_L) * (Omega_m/(ScaleFac*ScaleFac*ScaleFac) + Omega_L); 
		for(int ii = 0; ii < 3; ++ii)
			T[ii][ii] -= addoa;
#endif

		dsyevj3(T, Q, w);

		e0[dims.j*dims.k*i + dims.k*j + k] = w[0];
		e1[dims.j*dims.k*i + dims.k*j + k] = w[1];
		e2[dims.j*dims.k*i + dims.k*j + k] = w[2];
	}

#endif
	return;
}

#if 0
void do_minmax() {

	for(int i = 0; i < dims.i; ++i)
	for(int j = 0; j < dims.j; ++j)
	for(int k = 0; k < dims.k; ++k) {

		// iterate over the 26 neighboring cells
		psi_real phi0 = phi[ind3(i,j,k)];
		psi_real lmax = -1.0e99; 
		psi_real lmin = 1.0e99; 

		for(int ox = -1; ox <= 1; ++ox) 
		for(int oy = -1; oy <= 1; ++oy) 
		for(int oz = -1; oz <= 1; ++oz) {
			if(ox == 0 && oy == 0 && oz == 0) continue;
			int iflat = ind3(i+ox, j+oy, k+oz);
			if(phi[iflat] < lmin) lmin = phi[iflat];
			if(phi[iflat] > lmax) lmax = phi[iflat];
		}

		if(lmax > phi0 && lmin > phi0) {
			//printf("Found local minimum at %d %d %d, phi0 = %f, min = %f, max = %f\n",i,j,k,phi0,lmin,lmax);
			minima_phi.push_back(i);
			minima_phi.push_back(j);
			minima_phi.push_back(k);
		}
		else if(lmax < phi0 && lmin < phi0) {
			//printf("Found local maximum at %d %d %d, phi0 = %f, min = %f, max = %f\n",i,j,k,phi0,lmin,lmax);
			maxima_phi.push_back(i);
			maxima_phi.push_back(j);
			maxima_phi.push_back(k);
		}
	}

	for(int i = 0; i < dims.i; ++i)
	for(int j = 0; j < dims.j; ++j)
	for(int k = 0; k < dims.k; ++k) {

		// iterate over the 26 neighboring cells
		psi_real det0 = e0[ind3(i,j,k)]*e1[ind3(i,j,k)]*e2[ind3(i,j,k)];
		psi_real lmax = -1.0e99; 
		psi_real lmin = 1.0e99; 

		for(int ox = -1; ox <= 1; ++ox) 
		for(int oy = -1; oy <= 1; ++oy) 
		for(int oz = -1; oz <= 1; ++oz) {
			if(ox == 0 && oy == 0 && oz == 0) continue;
			int iflat = ind3(i+ox, j+oy, k+oz);
			psi_real dtmp = e0[iflat]*e1[iflat]*e2[iflat];
			if(dtmp < lmin) lmin = dtmp;
			if(dtmp > lmax) lmax = dtmp;
		}

		if(lmax > det0 && lmin > det0) {
			//printf("Found local minimum (det) at %d %d %d, phi0 = %f, min = %f, max = %f\n",i,j,k,det0,lmin,lmax);
			minima_det.push_back(i);
			minima_det.push_back(j);
			minima_det.push_back(k);
		}
		else if(lmax < det0 && lmin < det0) {
			//printf("Found local maximum (det) at %d %d %d, phi0 = %f, min = %f, max = %f\n",i,j,k,det0,lmin,lmax);
			maxima_det.push_back(i);
			maxima_det.push_back(j);
			maxima_det.push_back(k);
		}
	}


	return;
}

// read in a multidimensional hdf5 file
void read_hdf5(string filename, string fieldname) {
	vector<int> dims;
	HDFGetDatasetExtent(filename, fieldname, dims);
	dims.i = dims[0]; dims.j = dims[1]; dims.k = dims[2];
	HDFReadDataset(filename, fieldname, rho);
	return;
}

// read in a Fortran binary file
void read_fortran(string filename) {
	ifstream ifs;
	ifs.open(filename, ios::binary);
	if(!ifs.good()) {
		printf(" Failed to open file.\n");
		exit(0);
	}
	int blksize0;
	ifs.read((char*)&blksize0, sizeof(int));
	dims.i = (int)(pow((psi_real)blksize0/sizeof(float), 1.0/3.0) + 1.0e-3);
	dims.j = dims.i; dims.k = dims.i;

	// allocate a temp array
	float* tdata = (float*) malloc(blksize0);
	ifs.read((char*)tdata, blksize0);
	int blksize1;
	ifs.read((char*)&blksize1, sizeof(int));
	if(blksize0 != blksize1) {
		printf(" Bad file formatting.\n");
		exit(0);
	}
	ifs.close();

	// copy floats into the psi_real array
	rho.resize(dims.i*dims.j*dims.k);
	for(int i = 0; i < dims.i*dims.j*dims.k; ++i) {
		rho[i] = tdata[i];
	}
	free(tdata);

	return;
}

typedef struct {
	int npart[6];
	psi_real mass[6];
	psi_real time;
	psi_real redshift;
	int flag_sfr;
	int flag_feedback;
	unsigned int npartTotal[6];
	int flag_cooling;
	int num_files;
	psi_real BoxSize;
	psi_real Omega0;
	psi_real OmegaLambda;
	psi_real HubbleParam;
	int flag_stellarage;
	int flag_metals;
	unsigned int npartTotalHighWord[6];
	int flag_entropy_instead_u;
	int flag_psi_realprecision;
	char fill[56];
} Gadget2Header;

void read_gadget_header(string filename) {
	
	ifstream ifs;
	ifs.open(filename, ios::binary);
	if(!ifs.good()) {
		printf(" Failed to open file.\n");
		exit(0);
	}
	int blksize;
	Gadget2Header header;
	ifs.read((char*)&blksize, sizeof(int) );
	printf("blksize = %d\n", blksize);
	if(blksize != sizeof(Gadget2Header)) {
		printf(" Invalid Gadget snapshot format.\n");
		exit(0);
	}
	ifs.read((char*)&header, sizeof(Gadget2Header));
	ifs.close();

	// metadata
	//redshift = header.redshift;
	ScaleFac = header.time;
	BoxSize = header.BoxSize/1000.0; // need in Mpc/h
	HubbleParam = header.HubbleParam;
	Omega_m = header.Omega0;
	Omega_L = header.OmegaLambda;

	printf("   BoxSize = %lf (Mpc/h)\n", BoxSize);
	printf("   a = %lf\n", ScaleFac);
	printf("   H0 = %lf (100 km/(Mpc*s))\n", HubbleParam);
	printf("   O_m = %lf\n", Omega_m);
	printf("   O_L = %lf\n", Omega_L);
		
	return;
}

#endif



#endif
