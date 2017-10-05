from setuptools import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension(
    "psi", 
    sources = ["src/PSI.c", "src/mesh.c","src/grid.c", "src/geometry.c", "src/skymap.c",
        "src/cpsi.c", "src/refine.c", "src/rtree.c", "src/fft.c", "src/beamtrace.c"],
    extra_objects = ["include/mesh.h","include/grid.h", "include/geometry.h", "include/skymap.h",
        "include/refine.h", "include/rtree.h", "include/fft.h", "include/beamtrace.h"],
    define_macros=[('PYMODULE', None), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ('HAVE_FFTW', None)],
)

setup(
    name = 'the-phase-sheet-intersector-beta',
    version = '3.0.2',
    author = 'Devon Powell',
    author_email = 'devonmpowell1@gmail.com',
    ext_modules=[c_ext],
    include_dirs=['include']+numpy.distutils.misc_util.get_numpy_include_dirs(),
)
