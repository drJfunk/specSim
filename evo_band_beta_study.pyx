#Define the Band function
cimport cython

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)




@cython.cdivision(True)
cdef float  Band( float x, float A, float Ep, float alpha, float beta):

    cdef float cond = (alpha-beta)*Ep/(2+alpha)
    cdef float val
    if x < cond:

        val = A*( pow(x/100., alpha) * exp(-x*(2+alpha)/Ep) )
    else:
        val = A* ( pow( (alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return val




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
cdef float beta = -3.5




@cython.cdivision(True)
cpdef float evo(float ene, float t):
    
        #cdef float alpha = -1.
        #cdef float beta = -2.2
        
        cdef float Ep = 300.

        #cdef float maxFlux = 6.51795*1.
        #cdef float maxFlux = 6.51795*1.5444521
        #cdef float maxFlux = 6.51795*2.3853323
        #cdef float maxFlux = 6.51795*3.6840315
        #cdef float maxFlux = 6.51795*5.6898102
        #cdef float maxFlux = 6.51795*8.78763934
        #cdef float maxFlux = 6.51795*13.57208808
        #cdef float maxFlux = 6.51795*20.96144001
        #cdef float maxFlux = 6.51795*32.37394014
        cdef float maxFlux = 6.51795*50.
        
        
        cdef float A= maxFlux*1.

        
        cdef float val = Band(ene,A,Ep,alpha,beta)
    
        return val




