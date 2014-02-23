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
cpdef float PulseIntegrator(float t, float emin, float emax):

        cdef float val
        val = float(quad(evo,emin,emax,args=(t))[0])
        return val

        





###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import logspace

cdef float alpha = -1.
cdef float beta = -2.2

from peerModel.peerModel import peerModel


pm = peerModel("lum1.0gamma300.0tau1.0epl0.0ee0.1eb0.01ed0.5.txt")
pm.GenPhotonSpectrum(z=.01)






@cython.cdivision(True)
cpdef float evo(float ene, float t):
    
    ##This is a time-independent sim 

    
    
    
        
    cdef float val = pm.PhotonSpectrum(ene)
    
    return val




