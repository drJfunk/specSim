#Define the Band function
cimport cython

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)


cdef extern from "/Users/jburgess/Research/specSim/physPulse.h":
    float physPulse(float, float, float, float, float, float, float)







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





@cython.cdivision(True)
cpdef float evo(float ene, float t, p):


    cdef float etaT = .1
    cdef float etaR = 1.
    cdef float etaW = .1

    cdef float tf=2.
    cdef float u0=1.E5

    
    cdef float val = physPulse(t,ene,tf,u0,etaT,etaR,etaW)
    
    return val




