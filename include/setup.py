#!/usr/bin/env python

from distutils.core import setup, Extension
import os
import subprocess

os.environ['CXX'] = '/usr/bin/g++'
os.environ['CC'] = '/usr/bin/g++'
cwd = os.path.dirname(os.path.realpath(__file__))

ADMM_module = Extension('_ADMM',
			sources = [cwd+'/ADMM.i', cwd+'/ADMM.cpp',cwd+'/CVXcanon/src/CVXcanon.cpp',cwd+'/CVXcanon/src/BuildMatrix.cpp'
				,cwd+'/CVXcanon/src/EcosProblem.cpp',cwd+'/CVXcanon/src/LinOpOperations.cpp'],
				#,cwd+'/Solver/NetLasso.hpp',cwd+'/Solver/ModSquare.hpp'
				#,cwd+'/Solver/Square.hpp',cwd+'/Solver/ProximalMap.hpp'
				#,cwd+'/Solver/Edge.hpp',cwd+'/Solver/Node.hpp'
				#,cwd+'Solver/NodeVar.hpp'],
			include_dirs = [cwd+'/CVXCanon/include/Eigen/',cwd+'/ecos/include',cwd+'/ecos/external/SuiteSparse_config',cwd+'/Solver/'],
			library_dirs = ['ecos','Solver'],
			libraries = ['ecos','Solver'],
			swig_opts = ['-c++','-I '+cwd+'/ADMM.i -I '+cwd+'/CVXcanon/include/Eigen/src/Core/'],
			extra_compile_args = ['-Wall','-O2','-pg','-std=c++0x', '-fopenmp', '-ftree-vectorize'],
			extra_link_args = ['-lrt','-lgomp'] #required for older glibc versions
                           )

setup (name = 'ADMM',
       version = '0.1',
       author      = "Abhijit Sharang",
       description = """Simple swig example from docs""",
       ext_modules = [ADMM_module],
       py_modules = ["ADMM"],
       )
