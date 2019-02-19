'''

    Example 1: Voxelization, tetrahedral and curvilinear finite elements
    
    Voxelization (volume sampling) is controlled with 'grid', the integer dimensions of the target grid,
    and 'window', which sets the extent of the grid in coordinate-space.
    
    Change 'order' to 0, 1, or 2 to play with tetrahedral, trilinear, and triquadratic
    finite elements. 'tol' controls the accuracy of the voxelization in units of the voxel size. 
    Try varying it between 0 and 1.
    
    Fields to voxelize are requested via 'fields'. Setting 'm' to 'None' tells PSI to voxelize 
    the mass field to new numpy array.

'''

import PSI as psi
import helpers as hlp
import numpy as np


# =================================== #

# Try changing this to 0, 1, or 2!
order = 1

# Change these to play with the refinement behavior
reftol = 0.01
reflvlmax = 4

# =================================== #

# retrieve positions of a unit element
pos = np.array(hlp.unit_elements[order]).reshape((-1,3))
print ' - pos =', pos
vel = np.random.sample(pos.shape)
mass = np.ones(1)

# perturb the position randomly, keeping it within the unit box
pos -= 0.5
pos *= 0.5 + 0.4*np.random.sample(pos.shape)
pos += 0.5

# make a mesh, using the direct array input mode
# NOTE: connectivity must be 32-bit ints
conn = np.array([range(pos.shape[0])], dtype=np.int32)# TODO: connectivity array must be 2D 
pbox = (None, None) # TODO: this is still necessary to pass, even if the mesh is not periodic
print ' - Connectivity =', conn
mesh = psi.Mesh(loader='array', posar=pos, velar=vel, massar=mass, connar=conn, box=pbox)

print 'Connectivity done'

# create the Grid, specifying the resolution and projection window
ngrid = 3*(128,) 
win = ((0.,0.,0.),(1.,1.,1.))
grid = psi.Grid(type='cart', n=ngrid, window=win) 

print 'Grid done'

# call PSI.voxels()
psi.voxels(grid=grid, mesh=mesh, mode='density', refine_tolerance=reftol, refine_max_lvl=reflvlmax)

# print the error and show the figure
elemmass = np.sum(mesh.mass)
voxmass = np.sum(grid.fields["m"])
err = np.abs(1.0-voxmass/elemmass)
print 'Global error = %.10e' % err 
hlp.makeFigs(grid.fields['m'], title='Example 2: Voxelization of a cosmological density field')


