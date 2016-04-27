#!/usr/bin/env python

from distutils.core import setup, Extension
import os
import subprocess

os.environ['CXX'] = 'g++'
os.environ['CC'] = 'g++'
cwd = os.path.dirname(os.path.realpath(__file__))

ADMM_module = Extension('_ADMM',
			sources = [cwd+'/ADMM.i', cwd+'/ADMM.cpp',cwd+'/CVXcanon/src/CVXcanon.cpp',cwd+'/CVXcanon/src/BuildMatrix.cpp'
				,cwd+'/CVXcanon/src/EcosProblem.cpp',cwd+'/CVXcanon/src/LinOpOperations.cpp'],
			include_dirs = [cwd+'/CVXCanon/include/Eigen/',cwd+'/ecos/include',cwd+'/ecos/external/SuiteSparse_config'],
			library_dirs = ['ecos'],
			libraries = ['ecos'],
			swig_opts = ['-c++','-I '+cwd+'/ADMM.i'],
			extra_compile_args = ['-Wall','-O0','-pg','-std=c++0x', '-fopenmp'],
			extra_link_args = ['-lrt','-lgomp'] #required for older glibc versions
                           )

setup (name = 'ADMM',
       version = '0.1',
       author      = "Abhijit Sharang",
       description = """Simple swig example from docs""",
       ext_modules = [ADMM_module],
       py_modules = ["ADMM"],
       )
