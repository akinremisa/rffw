#%% *** This code is better runnning in Terminal***. '
# The code to recompile the fortran code in case any modifications 
# are made or required and then compiled using the f2py module in numpy
import numpy
from numpy import f2py

with open('rffw.f90') as sourcefile:
    sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='rffw', extension ='.f90', extra_args=['hrftn.f', 'sacsubf.f'])

#%% Run example code to test the compiled code.
from numpy import array
import matplotlib.pyplot as plt
import rffw
from numpy import arange

depth = array([ 5.0, 10.0, 20.0, 35.0, 50.0])
rho = array([2.67, 2.85, 2.95, 3.1, 3.3])
vp = array([4.0, 5.5, 7.0, 8.1, 8.4])
vs = array([2.8, 3.2, 3.6, 4.1, 4.5])
moddim = 5
rayp = 0.075
gaussalp = 2.0
delay = 5
n = 1024
delta = 0.04
rf = rffw.rffw(depth, vp, rho, vs, \
                rayp, gaussalp, delay, n, delta)
l = len(rf)
t = arange(0, l)
t = (delta *  t) - delay

plt.plot(t,rf)
plt.xlim(-delay, 40)
plt.show()
# %%