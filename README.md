# PSI

`PSI` (the Phase Space Intersector) is a Python module for post-processing the outputs of cosmological simulations. `PSI` applies the Phase Space Sheet discretization (see [Abel et al. 2012](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427...61A/abstract)) to a set of particles, then volume-samples the density and velocity fields (as well as higher-order moments) onto a Cartesian grid. Anyone can use `PSI` to produce smooth and continuous density and velocity maps from particle data!  The volume-sampling procedure is geometrically exact, and converves mass down to machine precision. 

The low-level volume-sampling routines are based on the [r3d software](https://github.com/devonmpowell/r3d),
described in 
[Powell & Abel (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563) and
[LA-UR-15-26964](https://github.com/devonmpowell/r3d/blob/master/la-ur-15-26964.pdf). 


---

## Usage:

- First, `git clone` the repository and `cd` into it 

- Then, `python setup.py install`. 

- To use in Python, simply `import PSI as psi` or similar 

- See `tutorial/Tutorial.ipynb` for a guide to using PSI! 

---

## License and attribution: 

PyPSI is free to use under the MIT License. See source files for the full license text.

If you use PyPSI for your research, please cite:

- Devon Powell and Tom Abel. An exact general remeshing scheme applied to physically conservative voxelization. Journal of Computational
Physics, 297:340–356, 2015a ([link](http://www.sciencedirect.com/science/article/pii/S0021999115003563))

- Devon M. Powell. r3d: Software for fast, robust geometric operations in 3D and 2D. Los Alamos National Lab technical report, LA‑UR‑15‑26964,
2015a ([link](https://github.com/devonmpowell/r3d/blob/master/la-ur-15-26964.pdf))


---

## Developers:

[Devon M. Powell](https://github.com/devonmpowell) (Stanford, MPA Garching)

[Tom Abel](https://github.com/yipihey) (Stanford, SLAC)

[Goran Jelic-Cizmek](https://github.com/JCGoran) (University of Geneva)

