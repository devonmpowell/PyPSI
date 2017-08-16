from distutils.core import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension("PSI", ["src/PSI.c", "src/mesh.c","src/grid.c", "src/geometry.c", "src/skymap.c",
    "src/cpsi.c", "src/refine.c", "src/rtree.c", "src/fft.c", "src/beamtrace.c"],
    define_macros=[('PYMODULE', None), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ('HAVE_FFTW', None)],
        )

setup(
    ext_modules=[c_ext],
    include_dirs=['include']+numpy.distutils.misc_util.get_numpy_include_dirs(),
)
