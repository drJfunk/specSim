#Define the Band function
cimport cython

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)







####Here is the pulse integrator
#### Uses GSL for fast integration

@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax,params):

        cdef float val
        val = float(quad(evo,emin,emax,args=(t,params))[0])
        return val







###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import logspace


from peerModel.peerModel import peerModel


pm = peerModel("/Users/jburgess/Research/specSim/lum0.1gamma100tau10epl0ee0.9eb1e-06ed0.1dtf500000.txt")
pm.GenPhotonSpectrum(z=1.)






@cython.cdivision(True)
cpdef float evo(float ene, float t, p):
    
    ##This is a time-independent sim 
        
    cdef float val = pm.PhotonSpectrum(ene)
    
    return val




