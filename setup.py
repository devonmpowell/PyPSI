from setuptools import setup, Extension
import numpy as np
import sys
import os


# default build behavior 
srcfiles = ["src/PSI.c", "src/mesh.c","src/grid.c", "src/geometry.c", "src/skymap.c",
        "src/cpsi.c", "src/refine.c", "src/rtree.c", "src/beamtrace.c",]
libs = []
libdirs = []
macros = [('PYMODULE', None), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),]


# determine the FFTW directories 
#fftw_home = "/usr/local/fftw3"
fftw_home = "/usr/local"
#fftw_home = "/scratch/users/dmpowel1/fftw-3.3.4"
print "FFTW3 home directory (Enter for default %s):"%fftw_home,
user_input = raw_input()
if not user_input.strip():
    user_input = fftw_home
if os.path.isfile(user_input+"/lib/libfftw3.a"):
    # if the fftw directory is good, add the path 
    fftw_home = user_input
    libs.append("fftw3")
    libdirs.append(fftw_home+"/lib")
    macros.append(('HAVE_FFTW', None))
    srcfiles.append("src/fft.c")
else:
    print "Invalid FFTW path %s (libfftw3.a not found)"%fftw_home
    user_input = ""
    while user_input not in ["Y", "n"]:
        print " Continue without building FFTW functionality? [Y/n]:",
        user_input = raw_input()
    if user_input == "n":
        print "Abort."
        sys.exit(0)


# define the C extension build 
c_ext = Extension(
    name = "PSI", 
    sources = srcfiles,
    define_macros = macros, 
    library_dirs = libdirs,
    libraries = libs,
)

# build the extension
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
