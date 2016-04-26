from distutils.core import setup
from subprocess import Popen
import shlex

#args = shlex.split('include/setup.py install')
#process = Popen(args)
#process.wait()

setup(name='snapvx',
    version='0.2',
    #data_files=[('Examples', [
    #    'Examples/BulkLoadData.csv',
    #    'Examples/BulkLoadEdges.edges',
    #    'Examples/BulkLoading.py',
    #    'Examples/HelloWorld.py',
    #    'Examples/LaplacianRegularization.py'])],
    py_modules=['snapvx']
    )


