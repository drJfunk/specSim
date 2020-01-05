
cimport cython

import numpy as np
cimport numpy as np

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)


cdef extern from "synchrotron.h":
    double synchrotron(double,double,double,double)
    double synchrotronComplex(double,double,double,double,double,double)
    double synchrotronPL(double, double, double, double, double)
    double synchrotronFast(double, double, double, double, double)
    double synchrotronPL_cutoff(double, double, double, double, double, double)
    double synchrotron_cutoff(double, double, double, double, double)
    double SSC(double, double, double , double )






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


Ep = 




@cython.cdivision(True)
cpdef float evo(float ene, float t):
    
               
        cdef float val
        val = synchrotron(ene,8,.01,3.5)
    
        return val




