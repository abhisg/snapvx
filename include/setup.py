#!/usr/bin/env python

from distutils.core import setup, Extension
import os

os.environ['CXX'] = 'g++-4.8'

ADMM_module = Extension('_ADMM',
			sources = ['ADMM_wrap.cxx', 'ADMM.cpp','CVXcanon/src/CVXcanon.cpp','CVXcanon/src/BuildMatrix.cpp'
				,'CVXcanon/src/EcosProblem.cpp','CVXcanon/src/LinOpOperations.cpp'],
			include_dirs = ['CVXCanon/include/Eigen/','ecos/include','ecos/external/SuiteSparse_config'],
			library_dirs = ['ecos'],
			libraries = ['ecos'],
                           )
"""ADMM_module = Extension('_ADMM',
			sources = ['ADMM_wrap.cxx', 'ADMM.cpp','CVXcanon/src/CVXcanon.cpp','CVXcanon/src/EcosProblem.cpp',
				'ecos/src/cone.c','ecos/src/ctrlc.c','ecos/src/ecos.c','ecos/src/equil.c','ecos/src/expcone.c',
				'ecos/src/kkt.c','ecos/src/preproc.c',
				'ecos/src/spla.c','ecos/src/splamm.c','ecos/src/timer.c','ecos/src/wright_omega.c'],
			include_dirs = ['CVXCanon/include','ecos/include','ecos/external/SuiteSparse_config','ecos/external/ldl/include/','ecos/external/amd/include/'],
			library_dirs = ['ecos'],
                           )"""

setup (name = 'ADMM',
       version = '0.1',
       author      = "Abhijit Sharang",
       description = """Simple swig example from docs""",
       ext_modules = [ADMM_module],
       py_modules = ["ADMM"],
       )
