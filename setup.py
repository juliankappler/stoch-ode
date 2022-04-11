import numpy, scipy
import os, sys
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True


if 'darwin' == (sys.platform).lower():
    extension = Extension("stoch_ode/*", ["stoch_ode/*.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-mmacosx-version-min=10.9'],
        extra_link_args=['-mmacosx-version-min=10.9'],
    )
else:
    extension = Extension("stoch_ode/*", ["stoch_ode/*.pyx"],
        include_dirs=[numpy.get_include()],
    )

setup(
    name='stoch_ODE',
    version='1.0.0',
    url='https://github.com/juliankappler/',
    author='Julian Kappler',
    license='MIT',
    description='python library for numerical simulation of Langevin equation',
    long_description='LangevinSim is a library for numerical simulation of the Langevin equation',
    platforms='works on all platforms (such as LINUX, macOS, and Microsoft Windows)',
    ext_modules=cythonize([ extension ],
        compiler_directives={"language_level": sys.version_info[0]},
        ),
    libraries=[],
    packages=['stoch_ode'],
    package_data={'stoch_ode': ['*.pxd']},
)
