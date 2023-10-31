#%% Setup the pip rffw package
import setuptools
with open("./README.md", "r") as fh:
    long_description = fh.read()
CLASSIFIERS = [
'Environment :: Console',
'Intended Audience :: Science/Research',
'License :: OSI Approved :: MIT License',
'Operating System :: OS Independent',
'Programming Language :: Python :: 3',
'Programming Language :: Python :: 3.7',
'Programming Language :: Python :: 3.8',
'Programming Language :: Python :: 3.9',
'Programming Language :: Python :: 3.10',
'Topic :: Scientific/Engineering :: Physics']

setuptools.setup(
    name="rffw",
    version="0.0.1",
    author='Stephen Akinremi, Islam Fadel',
    author_email=['s.akinremi@utwente.nl', 'i.e.a.m.fadel@utwente.nl'],
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    setup_requires=["numpy"],
    packages=setuptools.find_packages(),
    install_requires=['numpy'],
    classifiers=CLASSIFIERS)


###############################################################################################
# Compile the Fortran code
import numpy
from numpy import f2py
with open('rffw.f90') as sourcefile:
    sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='rffw', extension ='.f90', extra_args=['hrftn.f', 'sacsubf.f'])
