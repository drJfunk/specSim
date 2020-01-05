cimport cython
from math import floor
import numpy as np
#cimport numpy as np
from scipy.interpolate import interp1d


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
cdef float  Band( float x, float A, float Ep, float alpha, float beta):

    cdef float cond = (alpha-beta)*Ep
    cdef float val
    if x < cond:

        val = A*( pow(x/100., alpha) * exp(-x/Ep) )
    else:
        val = A* ( pow( (alpha -beta)*Ep/(100.),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return val        






eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t, p):



    cdef float val



    #if t> maxT or t<minT:
    #    return 0.
    


    val = Band(ene,p[0],p[1],p[2],p[3])

    
        
    return val




