import numpy as np
import matplotlib.pyplot as plt

# Makes a 2-panel plot with colorbar
def plot_twopanel_with_cb(dataslice, dataprojection, colsz=5.0, title=None, units=None, colors=plt.cm.RdBu_r):
    cbwidth = 0.03
    fig, (ax0, ax1, cax) = plt.subplots(1,3,figsize=((2.0+cbwidth)*colsz, colsz), gridspec_kw={'width_ratios': [1.0, 1.0, cbwidth]});
    plt.subplots_adjust(wspace=0.01)
    vmin = min(dataslice.min(), dataprojection.min())
    vmax = max(dataslice.max(), dataprojection.max())
    imargs = {'origin': 'lower', 'aspect': 'equal', 'interpolation': 'nearest', 'vmin': vmin, 'vmax': vmax, 'cmap': colors}
    im = ax0.imshow(dataslice, **imargs)
    ax0.set_xticks([]);
    ax0.set_yticks([]);
    ax0.set_title("Slice")
    im = ax1.imshow(dataprojection, **imargs)
    ax1.set_xticks([]);
    ax1.set_yticks([]);
    ax1.set_title("Projection")
    fig.colorbar(im, cax=cax)
    cax.set_ylabel(units)
    fig.suptitle(title, fontsize=14);
    return fig, (ax0,ax1), cax

# Makes a plot of positions and covariances in 2D 
def scatter_ellipses(mass, pos, cov=None,  title='', color='black', maxsz=25, ran=(0, 1, 0, 1), nstd=1):
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

# define sets of vertices 
unit_tetrahedron = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # tetrahedron (4 verts)
                             [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).reshape((-1,3))
unit_trilinear = np.array([[[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]], # trilinear (8 verts)
            [[[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]], [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]]]).reshape((-1,3))
unit_triquadratic = np.array([[[[0.0, 0.0, 0.0],[0.0, 0.0, 0.5], [0.0, 0.0, 1.0]], # triquadratic (27 verts) 
            [[0.0, 0.5 , 0.0], [0.0, 0.5 , 0.5], [0.0, 0.5 , 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.5], [0.0, 1.0, 1.0]]],
            [[[0.5,  0.0, 0.0], [0.5,  0.0, 0.5], [0.5,  0.0, 1.0]],
            [[0.5,  0.5 , 0.0], [0.5,  0.5 , 0.5], [0.5,  0.5 , 1.0]],
            [[0.5,  1.0, 0.0], [0.5,  1.0, 0.5], [0.5,  1.0, 1.0]]],
            [[[1.0, 0.0, 0.0], [1.0, 0.0, 0.5], [1.0, 0.0, 1.0]],
            [[1.0, 0.5 , 0.0], [1.0, 0.5 , 0.5], [1.0, 0.5 , 1.0]],
            [[1.0, 1.0, 0.0], [1.0, 1.0, 0.5], [1.0, 1.0, 1.0]]]]).reshape((-1,3))

# all unit elements together
unit_elements_pos = [unit_tetrahedron, unit_trilinear, unit_triquadratic]
unit_elements_conn = [np.arange(unit_tetrahedron.shape[0], dtype=np.int32).reshape(1,-1), 
                      np.arange(unit_trilinear.shape[0], dtype=np.int32).reshape(1,-1), 
                      np.arange(unit_triquadratic.shape[0], dtype=np.int32).reshape(1,-1)]

