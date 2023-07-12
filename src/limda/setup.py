from distutils.core import setup, Extension
from Cython.Build import cythonize
from numpy import get_include

ext = Extension("neighbor", sources=["neighbor.pyx"], include_dirs=['.', get_include()])
setup(name="neighbor", ext_modules=cythonize([ext]))
