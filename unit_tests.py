import unittest
from unittest import TestCase, main
import numpy as np
import PSI as psi

import sys

import matplotlib.pyplot as plt
import healpy as hp

# test basic PSIMOD functionality 
class PSIMODTests(TestCase):

    def test_phi_simple(self):

        return

        # loads, voxelizes, and computes the potential field
        print
        dim = 3
        box = 128 
        win = 21.0 
        G = 21.12131
        nk = 2
        print "g = ", G
        grid = psi.Grid(type='cart', n=dim*(box,), window=(dim*(0.,), dim*(win,))) 

        # set phi min, max such that phi min/max = -+1
        rho0 = -(np.pi*nk**2)/(G*win**2)
        grid.fields['m'][:,:,:] = rho0*np.sin(2*np.pi*nk*np.arange(box)/box)
        phi = psi.phi(grid, Gn=G)
        print 'phi min, max =', np.min(phi), np.max(phi)

        # plot the error
        phi -= np.sin(2*np.pi*nk*np.arange(box)/box)
        plt.imshow((phi[:,32,:]))
        plt.colorbar()
        plt.show()

    def test_phi_vox(self):

        return

        # loads, voxelizes, and computes the potential field
        print
        dim = 3
        box = 256 
        G = 12.120978
        mesh = psi.Mesh(filename='data/snapshot_999', loader='gadget2')
        #mesh = psi.Mesh(filename='data/snapshot_000', loader='gadget2')
        grid = psi.Grid(type='cart', n=dim*(box,), window=(mesh.boxmin, mesh.boxmax)) 
        psi.voxels(grid=grid, mesh=mesh)

        plt.imshow(np.log10(grid.fields['m'][:,32,:]))
        plt.colorbar()
        plt.show()


        phi = psi.phi(grid, Gn=G)
        plt.imshow((phi[:,32,:]))
        plt.colorbar()
        plt.show()

    def test_beam_tracing(self):


        #print
        grid = psi.Grid(type='hpring', n=32) 

        psi.beamtrace(grid=grid, metric='minkowski', obspos=(0.,0.,0.), obsvel=(0.,0.,0.))

        
        # show the pixel area plot
        hp.mollview(np.log10(grid.fields['m']), title='Mass map, err = %.5e'%err)
        plt.show()



    def test_annihilation(self):

        return

        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        #mesh = psi.Mesh(filename='', loader='hacky_test')
        #mesh = psi.Mesh(filename='data/box128_000', loader='gadget2')
        #print mesh.pos[mesh.connectivity][12475]

        grid = psi.Grid(type='hpring', n=32) 

        psi.skymap(grid=grid, mesh=mesh, bstep=2, mode=1)

        err = 1.0-np.sum(grid.fields["m"])
        print "mass = ", np.sum(grid.fields["m"]), 'err =', err

        #grid.fields["m"].tofile('data/scratch.np')
        
        # show the pixel area plot
        hp.mollview(np.log10(grid.fields['m']), title='Mass map, err = %.5e'%err)
        plt.show()

    def test_skymap(self):

        return

        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        #mesh = psi.Mesh(filename='', loader='hacky_test')
        #mesh = psi.Mesh(filename='data/box128_000', loader='gadget2')
        #print mesh.pos[mesh.connectivity][12475]

        grid = psi.Grid(type='hpring', n=32) 

        psi.skymap(grid=grid, mesh=mesh, bstep=2, mode=0)

        err = 1.0-np.sum(grid.fields["m"])
        print "mass = ", np.sum(grid.fields["m"]), 'err =', err

        #grid.fields["m"].tofile('data/scratch.np')
        
        # show the pixel area plot
        hp.mollview(np.log10(grid.fields['m']), title='Mass map, err = %.5e'%err)
        plt.show()

    def test_voxels(self):

        return

        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        #mesh = psi.Mesh(filename='data/box128_010', loader='gadget2')
        #print mesh.pos[mesh.connectivity][12475]

        grid = psi.Grid(type='cart', n=(64,64,64), window=(mesh.boxmin, mesh.boxmax)) 
        #grid = psi.Grid(type='cart', n=(64,64,64), window=((-5,-5,-5),(45,45,45))) 

        psi.voxels(grid=grid, mesh=mesh)

        #plt.imshow(np.log10(grid.fields["m"][:,:,32]))
        #plt.show()

        print "mass error = ", 1.0-np.sum(grid.fields["m"])
        
        # show the pixel area plot
        #err = grid.fields["m"]/grid.d-1.0
        #cran = max(np.abs(np.min(err)), np.abs(np.max(err)))
        #hp.mollview(err, title='Pix area error', cmap=plt.cm.RdYlGn,
                #min=-cran, max=cran)
        #plt.show()



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
        del grid

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
        del grid



# test basic Mesh functionality 
class MeshTests(TestCase):

    def test_mesh_init(self):
        print
        mesh = psi.Mesh(filename='data/snapshot_010', loader='gadget2')
        print ' - pos =', mesh.pos.dtype, mesh.pos.shape
        print ' - vel =', mesh.vel.dtype, mesh.vel.shape
        print ' - conn =', mesh.connectivity.dtype, mesh.connectivity.shape
        print ' - box =', mesh.boxmin, mesh.boxmax
        del mesh

if __name__ == '__main__':
    main(verbosity=2)
