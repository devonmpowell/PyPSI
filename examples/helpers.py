import numpy as np
import matplotlib.pyplot as plt

# Makes a 6-panel plot of projections and slices through the density field 
def makeFigs(voxdata, log=False, title='', colors=plt.cm.bone_r):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.subplots_adjust(wspace=0.0, hspace=0.0, bottom=0.08, left=0.05, top=0.92, right=0.85)
    axcb = fig.add_axes([0.86, 0.08, 0.03, 0.84])
    imargs = {'origin': 'lower', 'aspect': 'auto', 'interpolation': 'nearest', 'cmap': colors}
    if log:
        imargs['vmin'] = np.min(np.log10(voxdata)[voxdata > 0.0])
        imargs['vmax'] = np.max(np.log10(voxdata)[voxdata > 0.0])
    else:
        imargs['vmin'] = np.min(voxdata)
        imargs['vmax'] = np.max(voxdata)
    for ax, name in enumerate(['x', 'y', 'z']):
        if log:
            im = axes[0][ax].imshow(np.log10(voxdata.take(voxdata.shape[ax]/2, axis=ax)), **imargs)
            axes[1][ax].imshow(np.log10(np.mean(voxdata, axis=ax)), **imargs)
        else:
            im = axes[0][ax].imshow(voxdata.take(voxdata.shape[ax]/2, axis=ax), **imargs)
            axes[1][ax].imshow(np.mean(voxdata, axis=ax), **imargs)
        axes[0][ax].set_xticks([]);
        axes[0][ax].set_yticks([]);
        axes[1][ax].set_xticks([]);
        axes[1][ax].set_yticks([]);
        axes[1][ax].set_xlabel(name)
    fig.suptitle(title)
    axes[0][0].set_ylabel('Slice')
    axes[1][0].set_ylabel('Projection')
    fig.colorbar(im, cax=axcb)	
    plt.show()

# Makes a plot of positions and covariances in 2D 
def scatterEllipses(mass, pos, cov=None,  title='', color='black', maxsz=25, ran=(0, 1, 0, 1), nstd=1):
    fig, ax = plt.subplots(1,1, figsize=(10, 10))
    mmax = np.max(mass)
    if cov is None:
        ax.scatter(pos.T[0], pos.T[1], s=(maxsz*mass/mmax)**(2./3), lw=0, c=color, alpha=0.2)
    else:
        from matplotlib.patches import Ellipse
        for m, x, s in zip(mass, pos, cov):
            vals, vecs = np.linalg.eigh(s)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:,order]
            width, height = 2 * nstd * np.sqrt(vals) 
            theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
            ellip = Ellipse(xy=x, width=width, height=height, angle=theta, lw=0, fc=color, alpha=(m/mmax)**(1./3))
            ax.add_artist(ellip)
    ax.set_xlim([ran[0], ran[1]])
    ax.set_ylim([ran[2], ran[3]])
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig.suptitle(title)
    plt.show()

# randomly pick out some sample locations
def sampleSphere(r=1.0, center=[0.,0.,0.], nsamp=1):
    samps = np.random.randn(nsamp, 3)
    ilen = np.sum(samps**2, axis=1)**-0.5
    samps *= r*ilen[:,None]
    samps += np.array(center) 
    return samps

unit_elements = [
    np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), # tetrahedron (4 verts)
    np.array([[[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]], # trilinear (8 verts)
            [[[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]], [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]]]),
    np.array([[[[0.0, 0.0, 0.0],[0.0, 0.0, 0.5], [0.0, 0.0, 1.0]], # triquadratic (27 verts) 
            [[0.0, 0.5 , 0.0], [0.0, 0.5 , 0.5], [0.0, 0.5 , 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.5], [0.0, 1.0, 1.0]]],
            [[[0.5,  0.0, 0.0], [0.5,  0.0, 0.5], [0.5,  0.0, 1.0]],
            [[0.5,  0.5 , 0.0], [0.5,  0.5 , 0.5], [0.5,  0.5 , 1.0]],
            [[0.5,  1.0, 0.0], [0.5,  1.0, 0.5], [0.5,  1.0, 1.0]]],
            [[[1.0, 0.0, 0.0], [1.0, 0.0, 0.5], [1.0, 0.0, 1.0]],
            [[1.0, 0.5 , 0.0], [1.0, 0.5 , 0.5], [1.0, 0.5 , 1.0]],
            [[1.0, 1.0, 0.0], [1.0, 1.0, 0.5], [1.0, 1.0, 1.0]]]])
]

