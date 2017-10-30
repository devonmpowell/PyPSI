import unittest
from unittest import TestCase, main
import numpy as np
import PSI as psi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import healpy as hp
import sys

# set tolerances for floating-point checks
errtol = 1.0e-4


class MassMapTests(TestCase):

    def test_skymap(self):


        return

        print
        snap = 'data/snapshot_010'
        nside = 64 # 16 # 32
        mesh = psi.Mesh(filename=snap, loader='gadget2')
        #mesh = psi.Mesh(filename='', loader='hacky_test')
        grid = psi.Grid(type='hpring', n=nside) 
        psi.skymap(grid=grid, mesh=mesh, bstep=2, mode=0)

        # show the mass map and error
        elemmass = np.sum(mesh.mass)
        voxmass = np.sum(grid.fields["m"])
        err = np.abs(1.0-voxmass/elemmass)
        print ' - total mass =', voxmass 
        print " - mass error = ", err 
        hp.mollview(np.log10(grid.fields['m']), title='Mass map')
        plt.show()
        self.assertAlmostEqual(err, 0.0, delta=errtol)
 

    def test_voxels(self):

        #return

        # load a mesh, make a grid, and voxelize
        print
        ngrid = 3*(64,) 
        snap = 'data/snapshot_010'
        mesh = psi.Mesh(filename=snap, loader='gadget2')
        grid = psi.Grid(type='cart', n=ngrid, window=(mesh.boxmin, mesh.boxmax)) 
        psi.voxels(grid=grid, mesh=mesh)

        # check the total mass
        # show a picture
        elemmass = np.sum(mesh.mass)
        voxmass = np.sum(grid.fields["m"])
        err = np.abs(1.0-voxmass/elemmass)
        print ' - total mass =', voxmass 
        print " - mass error = ", err 
        plt.imshow(np.log10(grid.fields["m"][:,:,32]))
        plt.show()
        self.assertAlmostEqual(err, 0.0, delta=errtol)
       

# tests FFT stuff.. potential, power spectrum, etc 
class FFTTests(TestCase):

    def test_phi_simple(self):

        # loads, voxelizes, and computes the potential field
        print
        dim = 3
        box = 128 
        win = 21.0 
        G = 21.12131
        nk = 2
        grid = psi.Grid(type='cart', n=dim*(box,), window=(dim*(0.,), dim*(win,))) 

        # set phi min, max such that phi min/max = -+1
        rho0 = -(np.pi*nk**2)/(G*win**2)
        grid.fields['m'][:,:,:] = rho0*np.sin(2*np.pi*nk*np.arange(box)/box)
        phi = psi.phi(grid, Gn=G)
        analytic = np.sin(2*np.pi*nk*np.arange(box)/box)
        maxerr = np.max(np.abs(phi-analytic))
        print ' - max err in phi =', maxerr

        # test that the max error is small
        self.assertAlmostEqual(maxerr, 0.0, delta=errtol)


# test basic Grid functionality 
class GridTests(TestCase):


    def test_grid_cart(self):


        print
        grid = psi.Grid(type='cart', n=(10,10,10)) # default is 64^3 unit box
        print ' - type =', grid.type 
        print ' - fields =', grid.fields['m'].dtype, grid.fields['m'].shape
        print ' - window =', grid.winmin, grid.winmax
        print ' - n =', grid.n
        print ' - d =', grid.d
        center, boundary, vol = grid.getCellGeometry(cell=[123])
        #print center.shape, boundary.shape, vol.shape
        #print center, boundary, vol

    def test_grid_hpring(self):
        print
        grid = psi.Grid(type='hpring')
        print ' - type =', grid.type 
        print ' - fields =', grid.fields['m'].dtype, grid.fields['m'].shape
        print ' - window =', grid.winmin, grid.winmax
        print ' - n =', grid.n
        print ' - d =', grid.d

        # get a pixel boundary
        center, boundary, vol = grid.getCellGeometry(cell=None, bstep=2)
        #print center, boundary, vol

# test basic Mesh functionality 
class MeshTests(TestCase):

    def test_mesh_init(self):
        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        print ' - pos =', mesh.pos.dtype, mesh.pos.shape
        print ' - vel =', mesh.vel.dtype, mesh.vel.shape
        print ' - conn =', mesh.connectivity.dtype, mesh.connectivity.shape
        print ' - box =', mesh.boxmin, mesh.boxmax

# test basic PSIMOD functionality 
class PSIMODTests(TestCase):

    def test_beam_tracing(self):

        #return



        #metric = psi.Metric(type='flrw', filepattern='/home/devon/HDD/PS-128-vm/geo/%04d')
        metric = psi.Metric(type='kerr')


        with open('/home/devon/HDD/PS-128-vm/geo/0137', 'rb') as f: 

            # read potential
            f.seek(0)
            phi = np.fromfile(f, count=128**3,
                dtype=np.float32).reshape((128,128,128))
            print phi.shape, np.min(phi), np.max(phi)

            # read velocity 
            f.seek(128**3)
            vel = np.fromfile(f, count=3*128**3,
                dtype=np.float32).reshape((128,128,128,3))

            # read potential
            f.seek(4*128**3)
            gradphi = np.fromfile(f, count=3*128**3,
                dtype=np.float32).reshape((128,128,128,3))
            print gradphi.shape, np.min(gradphi), np.max(gradphi)

            f.close()

        gradax = 2
        mygrad = np.roll(phi,1,axis=gradax)-np.roll(phi,-1,axis=gradax)

        fig, ax = plt.subplots(1,2,figsize=(20,10))
        #ax[0].imshow(phi[64,:,:])


        ax[0].imshow(mygrad[64,:,:])
        ax[1].imshow(gradphi[64,:,:,gradax])
        #ax[1].contour(phi[64,:,:])
        plt.show()


        #return

        #print
        grid = psi.Grid(type='hpring', n=16) 
        
        rayinfo = psi.beamtrace(grid=grid, metric=metric, obspos=(10.,9.,8.), obsvel=(0.,0.,0.))

        print rayinfo.shape
        print rayinfo[0,:,:]



        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for tray in rayinfo:
            ray = tray[tray[:,0] > 0.0]
            ax.plot(ray[:,1], ray[:,2], ray[:,3], 'k-', lw=0.2, label='parametric curve')
        #ax.set_xlim(-30,30)
        #ax.set_ylim(-30,30)
        #ax.set_zlim(-30,30)
        plt.show()



        # show the pixel area plot
        hp.mollview((grid.fields['m']), title='minkowski')
        plt.show()



    def test_annihilation(self):

        return

        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        #mesh = psi.Mesh(filename='', loader='hacky_test')
        #mesh = psi.Mesh(filename='data/box128_010', loader='gadget2')
        #print mesh.pos[mesh.connectivity][12475]

        grid = psi.Grid(type='hpring', n=128) 

        psi.skymap(grid=grid, mesh=mesh, bstep=1, mode=1)

        err = 1.0-np.sum(grid.fields["m"])
        print "mass = ", np.sum(grid.fields["m"]), 'err =', err
        print " - min = ", np.min(grid.fields["m"])
        print " - max = ", np.max(grid.fields["m"])

        grid.fields["m"].tofile('data/annihilation128_shell.np')
        
        # show the pixel area plot
        #hp.mollview(np.log10(grid.fields['m']), title='Mass map, err = %.5e'%err)
        hp.mollview(np.log10(grid.fields['m']))#, title='Mass map, err = %.5e'%err)#, min = -4)
        plt.show()



####### runs all tests by default #######
if __name__ == '__main__':
    main(verbosity=2)
