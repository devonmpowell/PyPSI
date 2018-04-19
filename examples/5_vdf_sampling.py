'''

    Example 5: VDF sampling  

    Demonstrates how to use PSI.VDF() and PSI.crossStreams()

'''
        
import PSI as psi
import helpers as hlp
import numpy as np


# load the mesh 
snap = '../data/snapshot_010'
mesh = psi.Mesh(filename=snap, loader='gadget2')

# make an RStarTree out of the mesh
rst = psi.RStarTree(mesh)

# define the velocity-dependent rate coefficient
def xsn(v):
    return 1.0 

# loop over all sample locations
for sample_pos in hlp.sampleSphere(r=0.2, center=[10.,10.,13.], nsamp=5): 

    print 'Next sample...'

    # sample the multi-streaming region
    rho, vel = psi.VDF(mesh, tuple(sample_pos), tree=rst)

    # compute bulk quantities from the stream data
    # the annihilation rate is computed from from xsnfunc()
    nstr, rhotot, vbulk, sigv, rate = psi.crossStreams(rho, vel, xsnfunc=xsn)

    # print out the computed data
    print ' - nstreams =', nstr
    print ' - rhotot =', rhotot
    print ' - vbulk =', vbulk
    print ' - dispersion =', sigv
    print ' - annihilation rate =', rate 
    print
