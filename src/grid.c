/*
 *
 * grid.c 
 *
 * Copyright (c) 2013-2023 Devon M. Powell
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "grid.h"
#include "geometry.h"

// declarations internal to this file
void neighbors(int order_, int pix, int* result);
void xyf2loc(psi_real x, psi_real y, psi_int face, psi_real* vecout);
int xyf2ring(int order_, int ix, int iy, int face_num);
void ring2xyf(int order_, int pix, int* ix, int* iy, int* face_num);

psi_int psi_grid_get_cell_geometry(psi_grid* grid, psi_dvec grind, psi_int step, psi_rvec* boundary, psi_rvec* center, psi_real* vol) {

	psi_int ii, jj, ix, iy, i, j, k, m, face;
	psi_int order_, nside_, npix_;
	psi_real dx[2], x0[2], xc, yc, dc, d;

	switch(grid->type) {
	
		case PSI_GRID_CART:

			printf("cart3d vertices for cell %d %d %d\n", grind.i, grind.j, grind.k);
		
			// back out the grid index
			// compute cellel boundary and center
			//ii = cell/grbp->idim[1]; // row-major order
			//jj = cell%grbp->idim[1]; // switch the / and % to change to column-major
			//dx[0] = grbp->fdim[0]/grbp->idim[0];
			//dx[1] = grbp->fdim[1]/grbp->idim[1];
			//x0[0] = -0.5*grbp->fdim[0];
			//x0[1] = -0.5*grbp->fdim[1];
			//memset(boundary, 0, 12*sizeof(real));
			//memset(center, 0, 3*sizeof(real));
	
			//// TODO: z-coordinates
			//boundary[3*0+0] = x0[0]+(ii+0)*dx[0];
			//boundary[3*0+1] = x0[1]+(jj+0)*dx[1];
			//boundary[3*1+0] = x0[0]+(ii+1)*dx[0];
			//boundary[3*1+1] = x0[1]+(jj+0)*dx[1];
			//boundary[3*2+0] = x0[0]+(ii+1)*dx[0];
			//boundary[3*2+1] = x0[1]+(jj+1)*dx[1];
			//boundary[3*3+0] = x0[0]+(ii+0)*dx[0];
			//boundary[3*3+1] = x0[1]+(jj+1)*dx[1];
	
			// average the center vertex for convenience
			//center[0] = 0.25*(boundary[3*0+0]+boundary[3*1+0]+boundary[3*2+0]+boundary[3*3+0]);
			//center[1] = 0.25*(boundary[3*0+1]+boundary[3*1+1]+boundary[3*2+1]+boundary[3*3+1]);
			//center[2] = 0.25*(boundary[3*0+2]+boundary[3*1+2]+boundary[3*2+2]+boundary[3*3+2]);
			break;

		case PSI_GRID_HPRING:

	  		order_  = grid->n.i; 
	  		nside_  = grid->n.j;
	  		npix_  = grid->n.k;
			if(grind.i >= npix_) return 0;
			ring2xyf(order_, grind.i, &ix, &iy, &face); // only RING ordering for now
			dc = 0.5 / nside_;
			xc = (ix + 0.5)/nside_, yc = (iy + 0.5)/nside_;
			d = 1.0/(step*nside_);
			xyf2loc(xc, yc, face, (psi_real*)center);
			for(i=0, j=step, k=2*step, m=3*step; i < step; ++i, ++j, ++k, ++m) {
			    xyf2loc(xc+dc-i*d, yc+dc, face, (psi_real*)&boundary[i]);
			    xyf2loc(xc-dc, yc+dc-i*d, face, (psi_real*)&boundary[j]);
			    xyf2loc(xc-dc+i*d, yc-dc, face, (psi_real*)&boundary[k]);
			    xyf2loc(xc+dc, yc-dc+i*d, face, (psi_real*)&boundary[m]);
			}
			*vol = 0.0;
			for(i = 0; i < step*4; ++i) {
				*vol += psi_omega3(boundary[i], boundary[(i+1)%(step*4)], *center);
			}

			break;


		default:
			return 0;
	}




	return 1;
}





//////////////////////////////////////////////////////
//       Healpix functions ported from the C++ beta
// 		https://github.com/healpy/healpy/tree/afe5a31fff01ff3680e62d550934ad7efca73040/hpbeta
//////////////////////////////////////////////////////

#define isqrt(x) ((int)sqrt(x+0.5)) // index matrices
const static int jrll[] = {2,2,2,2,3,3,3,3,4,4,4,4};
const static int jpll[] = {1,3,5,7,0,2,4,6,1,3,5,7};
const static int xoffset[] = { -1,-1, 0, 1, 1, 1, 0,-1 };
const static int yoffset[] = {  0, 1, 1, 1, 0,-1,-1,-1 };
const static int facearray[][12] =
      { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },   // S
        {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },   // SE
        { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },   // E
        {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },   // SW
        {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },   // center
        {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },   // NE
        { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },   // W
        {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },   // NW
        {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } }; // N
const static int swaparray[][12] =
      { {  0,0,0,0,0,0,0,0,3,3,3,3 },   // S
        {  0,0,0,0,0,0,0,0,6,6,6,6 },   // SE
        {  0,0,0,0,0,0,0,0,0,0,0,0 },   // E
        {  0,0,0,0,0,0,0,0,5,5,5,5 },   // SW
        {  0,0,0,0,0,0,0,0,0,0,0,0 },   // center
        {  5,5,5,5,0,0,0,0,0,0,0,0 },   // NE
        {  0,0,0,0,0,0,0,0,0,0,0,0 },   // W
        {  6,6,6,6,0,0,0,0,0,0,0,0 },   // NW
        {  3,3,3,3,0,0,0,0,0,0,0,0 } }; // N

void neighbors(int order_, int pix, int* result) {

  // back out all the Healpix parameters
  int nside_  = 1<<order_;
int ix, iy, face_num, m, i;
//(scheme_==RING) ?
  ring2xyf(order_, pix,&ix,&iy,&face_num);// : nest2xyf(pix,ix,iy,face_num);

const int nsm1 = nside_-1;
if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
  {
  //if (scheme_==RING)
    for (m=0; m<8; ++m)
      result[m] = xyf2ring(order_, ix+xoffset[m],iy+yoffset[m],face_num);
  //else
    //for (int m=0; m<8; ++m)
      //result[m] = xyf2nest(ix+xoffset[m],iy+yoffset[m],face_num);
  }
else
  {
  for (i=0; i<8; ++i)
    {
    int x=ix+xoffset[i];
    int y=iy+yoffset[i];
    int nbnum=4;
    if (x<0)
      { x+=nside_; nbnum-=1; }
    else if (x>=nside_)
      { x-=nside_; nbnum+=1; }
    if (y<0)
      { y+=nside_; nbnum-=3; }
    else if (y>=nside_)
      { y-=nside_; nbnum+=3; }

    int f = facearray[nbnum][face_num];
    if (f>=0)
      {
      if (swaparray[nbnum][face_num]&1) x=nside_-x-1;
      if (swaparray[nbnum][face_num]&2) y=nside_-y-1;
      if (swaparray[nbnum][face_num]&4) {
	  	int tmp = x;
		x = y;
		y = tmp;
	  } 
      //result[i] = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
      result[i] = xyf2ring(order_, x,y,f);
      }
    else
      result[i] = -1;
    }
  }
}


// returns the normalized Cartesian coordinate given an x,y coordinate for the given face
void xyf2loc(psi_real x, psi_real y, psi_int face, psi_real* vecout) {
	psi_real z, phi, sinTheta, tmp;
    psi_real jr = jrll[face] - x - y;
    psi_real nr;
    if (jr<1) {
        nr = jr;
        tmp = nr*nr/3.;
        z = 1 - tmp;
    } else if (jr>3) {
      nr = 4-jr;
      tmp = nr*nr/3.;
      z = tmp - 1;
    } else {
        nr = 1;
        z = (2-jr)*2./3.;
    }
    tmp=jpll[face]*nr+x-y;
    if (tmp<0) tmp+=8;
    if (tmp>=8) tmp-=8;
    phi = (nr<1e-15) ? 0 : (0.25*PI*tmp)/nr;

	// convert to Cartesian
	sinTheta = sqrt((1.0-z)*(1.0+z));
	vecout[0] = sinTheta*cos(phi);
	vecout[1] = sinTheta*sin(phi);
	vecout[2] = z;
}

int xyf2ring(int order_, int ix, int iy, int face_num) {

  int nside_  = 1<<order_;
  int npface_ = nside_<<order_;
  int ncap_   = (npface_-nside_)<<1;
  int npix_   = 12*npface_;

  int nl4 = 4*nside_;
  int jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  int nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = npix_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }


// gets the face and pixel indices from one NEST index
void ring2xyf(int order_, int pix, int* ix, int* iy, int* face_num) {
  int iring, iphi, kshift, nr;

  // back out all the Healpix parameters
  int nside_  = 1<<order_;
  int nl2 = 2*nside_;
  int npface_ = nside_<<order_;
  int ncap_   = (npface_-nside_)<<1;
  int npix_   = 12*npface_;

  if (pix<ncap_) // North Polar cap
    {
    iring = (1+isqrt(1+2*pix))>>1; //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    *face_num=(iphi-1)/nr;
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    int ip = pix - ncap_;
    int tmp = (order_>=0) ? ip>>(order_+2) : ip/(4*nside_);
    iring = tmp+nside_;
    iphi = ip-tmp*4*nside_ + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
   	int ire = iring-nside_+1,
      irm = nl2+2-ire;
    int ifm = iphi - ire/2 + nside_ -1,
      ifp = iphi - irm/2 + nside_ -1;
    if (order_>=0)
      { ifm >>= order_; ifp >>= order_; }
    else
      { ifm /= nside_; ifp /= nside_; }
    *face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else // South Polar cap
    {
    int ip = npix_ - pix;
    iring = (1+isqrt(2*ip-1))>>1; //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    *face_num = 8 + (iphi-1)/nr;
    }

  int irt = iring - (jrll[*face_num]*nside_) + 0;
  int ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy = (-ipt-irt) >>1;
}

