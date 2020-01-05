cimport cython
from math import floor
import numpy as np
#cimport numpy as np
from scipy.interpolate import interp1d


from evolver import Evolver

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)


@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax, params):
    cdef float val
    val = float(quad(evo,emin,emax,args=(t,params))[0])
    return val

        

@cython.cdivision(True)
cdef float  PL( float x, float A, float index):

    
    return A*pow(x/100.,index)




#minT = tStart.min()
#maxT = tStop.max()



eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax
evol = Evolver()

@cython.cdivision(True)
cpdef float evo(float ene, float t, p):



    cdef float val




    evol.SetPulse(p[0],p[1],p[2],p[3])
    evol.SetSpectralEvolution(p[4],p[5],p[6])
    evol.SetSpectrum(p[7],p[8])

    val = evol.IntegratedSpectrum(ene,p[9],p[10])

    
        
    return val




