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
cdef float  CPL( float x, float A, float index, float Ep):


    cdef float val

    val = A*pow(x/300.,index)*exp(-x/Ep)

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





@cython.cdivision(True)
cpdef float evo(float ene, float t):
    
        #cdef float alpha = -1.
        #cdef float beta = -2.2
        
        cdef float Ep = 300.

        
        #cdef float maxFlux = 2.17867331*1.
        #cdef float maxFlux = 2.17867331*1.5444521
        #cdef float maxFlux = 2.17867331*2.3853323
        #cdef float maxFlux = 2.17867331*3.6840315
        #cdef float maxFlux = 2.17867331*5.6898102
        #cdef float maxFlux = 2.17867331*8.78763934
        #cdef float maxFlux = 2.17867331*13.57208808
        #cdef float maxFlux = 2.17867331*20.96144001
        #cdef float maxFlux = 2.17867331*32.37394014
        cdef float maxFlux = 2.17867331*50.
        
        
        cdef float A= maxFlux*1.

        
        cdef float val = CPL(ene,A,alpha,Ep)
    
        return val




