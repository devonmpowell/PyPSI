'''
    
    Example 4: RStarTree 

    Demonstrates how to use the RStarTree class
    Still somewhat rudimentary

    RStarTree is used by PSI.VDF and PSI.voxels(mode='annihlation') so far.

'''
        
import PSI as psi
import helpers as hlp
import numpy as np


# load the mesh 
snap = '../data/snapshot_010'
mesh = psi.Mesh(filename=snap, loader='gadget2')

# make an RStarTree out of the mesh
# Currently this is the only constructor;
# a mesh must be provided
rst = psi.RStarTree(mesh)

# set up the query box
query_box = ( (10., 9., 8.), (10.3, 9.5, 8.2) )

# query the tree
# this returns an index for each queried element!
qinds = rst.query(box=query_box)
print ' - element indices =', qinds

# use fancy indexing to get the element connectivity back
print ' - vertex ids =', mesh.connectivity[qinds]
print ' - element coordinates =', mesh.pos[mesh.connectivity[qinds]]
