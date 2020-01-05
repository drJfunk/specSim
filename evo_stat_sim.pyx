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
def pl(x,A,index):
    return A * (x/300.)**(index)



eMin = 7.5
eMax = 40000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t,p):


      val = pl(ene,p[0],-2.)
      return val




