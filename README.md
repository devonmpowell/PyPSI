## PSI: The Phase Space Intersector 

A Python extension for volume-sampling and point-sampling cosmological density and velocity fields.

The volume-sampling kernel is based on the [r3d software](https://github.com/devonmpowell/r3d),
described in 
[Powell & Abel (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563) and
[LA-UR-15-26964](la-ur-15-26964.pdf). 

---

### Features:

- PSI! 

- Example code in `examples/`

---

### Usage:

- First, `git clone` the repository and `cd` into it 

- Then, `python setup.py install`. You will be prompted for the location of the [FFTW3 library](http://www.fftw.org/); 
if you do not have a copy this is OK. PyPSI will simply be compiled without functions that use
FFTs.

- To use in Python, simply `import PSI as psi` or similar 


---

### Licensing: 

Coming soon.

For now, use as you wish for non-commercial purposes with proper attribution.


---

### Developers:

Devon Powell (Stanford, MPA Garching)

Tom Abel (Stanford, SLAC)

Goran Jelic-Cizmek (University of Geneva)
