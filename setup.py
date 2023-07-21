from setuptools import setup, Extension
import numpy as np
import sys
import os


# default build behavior
srcfiles = ["src/PSI.c", "src/mesh.c", "src/grid.c", "src/geometry.c",
        "src/cpsi.c", "src/refine.c", "src/rtree.c",]
libs = []
libdirs = []
macros = [('PYMODULE', None), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),]


# define the C extension build
c_ext = Extension(
    name = "PSI",
    sources = srcfiles,
    define_macros = macros,
    library_dirs = libdirs,
    libraries = libs,
    extra_compile_args = ["-w"],
)

# build the extension
setup(
    name = 'the-phase-space-intersector',
    version = '3.1.2',
    author = 'Devon Powell',
    author_email = 'devonmpowell1@gmail.com',
    url = 'https://github.com/devonmpowell/PyPSI',
    ext_modules=[c_ext],
    install_requires = ['numpy'],
    include_dirs=['include', np.get_include()],
)
