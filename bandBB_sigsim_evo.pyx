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

    cdef float cond = (alpha-beta)*Ep/(2+alpha)
    cdef float val
    if x < cond:

        val = ( pow(x/100., alpha) * exp(-x*(2+alpha)/Ep) )
    else:
        val =  ( pow( (alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return A*val





def bb(x,A,kT):
    return A * x**2 * (exp(x/kT)-1.)**-1.



eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t,p):



    val = 5.*( Band(ene,0.01,500.,-1.,-2.2)+bb(ene,1.E-6,40.) )
    
    
        
    return val




