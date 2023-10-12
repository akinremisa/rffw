"""
rffw is a python module for forward calculation of 1D radial receiver functions. 
It is based on the hrftn96 program (Computer programs in seismology; Herrmann, 2013). 
Th original program is written in Fortran and is available at http://www.eas.slu.edu/eqc/eqccps.html.
rffw is the python module that compiles the Fotran code into a python lirary using f2py.

# Installation
Here are some instructions to install rf into a fresh conda environment::

1. Deactivate running Conda Env
conda deactivate

2. Create new conda environment 'rffw' 
conda create -n rffw python=3.9 anaconda -y

3. Activate the environement 
conda activate rffw

4. Download via wget or web-browser
wget https://github.com/akinremisa/rffw/archive/refs/heads/main.zip

5. unzip the package
unzip main.zip  

6. Change directory to the same directory that this repo is in (i.e., same directory as setup.py)
cd rffw-main 
 
7. run pip comand:
pip install - e. 


Alternatively, the environment can be created using the environment.yml file in the root directory of this repo::
conda env create -f environment.yml



References
Herrmann, R. B. (2013). Computer programs in seismology: An evolving tool for
instruction and research. Seismological Research Letters, 84 , 1081-1088. doi: 10.1785/0220110096
"""