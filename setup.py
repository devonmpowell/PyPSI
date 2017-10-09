from setuptools import setup, Extension
import numpy as np

c_ext = Extension(
    name = "PSI", 
    sources = ["src/PSI.c", "src/mesh.c","src/grid.c", "src/geometry.c", "src/skymap.c",
        "src/cpsi.c", "src/refine.c", "src/rtree.c", "src/fft.c", "src/beamtrace.c"],
    extra_objects = ["include/mesh.h","include/grid.h", "include/geometry.h", "include/skymap.h",
        "include/refine.h", "include/rtree.h", "include/fft.h", "include/beamtrace.h"],
    define_macros = [('PYMODULE', None), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ('HAVE_FFTW', None)],
)

setup(
    name = 'the-phase-space-intersector',
    version = '3.0.2',
    author = 'Devon Powell',
    author_email = 'devonmpowell1@gmail.com',
    url = 'https://upload.pypi.org/legacy/',
    ext_modules=[c_ext],
    install_requires = ['numpy'],
    include_dirs=['include', np.get_include()],
)
