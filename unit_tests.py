import unittest
from unittest import TestCase, main
import numpy as np
import PSI as psi


import matplotlib.pyplot as plt
import healpy as hp

# test basic PSIMOD functionality 
class PSIMODTests(TestCase):

    def test_skymap(self):

        print
        mesh = psi.Mesh(filename='data/snapshot_010', type='gadget2')
        grid = psi.Grid(type='hpring', n=16) 

        psi.skymap(grid=grid, mesh=mesh, bstep=1)
        print "skymap n=", grid.n

        hp.mollview(grid.fields["m"]/grid.d, title='Pix area fraction')
        plt.show()
        #print grid.fields["m"]
        #print grid.fields["m"]


# test basic Grid functionality 
class GridTests(TestCase):

    def test_grid_cart3d(self):
        print
        grid = psi.Grid(type='cart3d', n=(10,10,10)) # default is 64^3 unit box
        print ' - type =', grid.type 
        print ' - fields =', grid.fields['m'].dtype, grid.fields['m'].shape
        print ' - window =', grid.winmin, grid.winmax
        print ' - n =', grid.n
        print ' - d =', grid.d
        #print grid.getCellGeometry(cell=[123])
        center, boundary, vol = grid.getCellGeometry(cell=[123])
        print center.shape, boundary.shape, vol.shape
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
        center, boundary, vol = grid.getCellGeometry(cell=[123, 456], bstep=2)
        print center
        print boundary
        print vol
        del grid



# test basic Mesh functionality 
class MeshTests(TestCase):

    def testMeshInit(self):
        print
        mesh = psi.Mesh(filename='data/snapshot_010', type='gadget2')
        print ' - type =', mesh.type
        print ' - file =', mesh.filename
        print ' - pos =', mesh.pos.dtype, mesh.pos.shape
        print ' - vel =', mesh.vel.dtype, mesh.vel.shape
        del mesh

if __name__ == '__main__':
    main(verbosity=2)
