## PSI: The Phase Space Intersector 

A Python extension for volume-sampling and point-sampling cosmological density and velocity fields.

The volume-sampling kernel is based on the [r3d software](https://github.com/devonmpowell/r3d),
described in 
[Powell & Abel (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563) and
[LA-UR-15-26964](la-ur-15-26964.pdf). 

---

### Features:

- PSI! 

---

### Usage:

- First, `git clone` the repository and `cd` into it 

- Then, `python setup.py install`. You will be prompted for the location of the [FFTW3 library](http://www.fftw.org/); 
if you do not have it installed this is OK. PyPSI will be compiled without functions that use
FFTs.

- To use in Python, simple `import PyPSI` 


---

### Licensing: 

???
