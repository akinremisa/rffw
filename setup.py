#%% Setup the pip rffw package
from numpy import f2py
with open('rffw.f90') as sourcefile:
    sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='rffw', extension ='.f90', extra_args=['hrftn.f', 'sacsubf.f'])


import setuptools

with open("./README.md", "r") as fh:
    long_description = fh.read()

CLASSIFIERS = ['Environment :: Console',
'Intended Audience :: Science/Research',
'License :: OSI Approved :: MIT License',
'Operating System :: OS Independent',
'Programming Language :: Python :: 3',
'Programming Language :: Python :: 3.7',
'Programming Language :: Python :: 3.8',
'Programming Language :: Python :: 3.9',
'Programming Language :: Python :: 3.10',
'Topic :: Scientific/Engineering :: Physics']
REQUIRES=['numpy<=1.23.5', 'cython>=0.29.21', 'blosc2~=2.0.0',
          'FuzzyTM>=0.4.0', 'OWSLib', 'geographiclib', 'pyshp>=2.1',
          'pyproj>=3.0.0']
setuptools.setup(
    name="rffw",
    version="0.0.1",
    author='Stephen Akinremi, Islam Fadel',
    author_email=['s.akinremi@utwente.nl', 'i.e.a.m.fadel@utwente.nl'],
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=REQUIRES,
    classifiers=CLASSIFIERS,
    python_requires='>=3.6')
