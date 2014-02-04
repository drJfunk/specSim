cimport cython
from numpy import linspace
from evo cimport evo, PulseIntegrator

from libc.stdlib cimport rand, RAND_MAX
cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)



#cdef class Function:
#    cpdef evaluate(self, float x):
#        return 0

  




@cython.cdivision(True)
def Sample(float t, float fMax, float xMin, float xMax):

        cdef int flag = 1
        cdef float energyGuess
        cdef float test
        cdef float fCheck
        #cdef float randnum
        while flag==1:
           
            energyGuess = xMin +  float(rand()) / float(RAND_MAX/(xMax - xMin))
            test = float(rand()) / float(RAND_MAX/(fMax))
            fCheck = evo(energyGuess,t)

            if  test <= fCheck:
                    flag = 0
        
        return energyGuess

@cython.cdivision(True)
def homoGen(float t0, float tMax,  float rate):
    cdef float t
    cdef float lastT
    times=[t0]
    while times[-1]<tMax:
        lastT = times[-1]
        t = lastT-(1./rate)*log( float(rand()) / float(RAND_MAX)  )
        times.append(t)
        
    return times



@cython.cdivision(True)
@cython.boundscheck(False)
def nonHomoGen(float t0, float tMax, float emin, float emax, float areaFrac):
        '''
        Non-homogeneous poisson process generator
        for a given max rate and time range, this function 
        generates time tags sampled from the energy integrated
        lightcurve.

        '''
        tmpT=linspace(t0,tMax,1000.)
        tmp = map(lambda x: PulseIntegrator(x,emin,emax)*areaFrac,tmpT)
        fmax = max(tmp)



        cdef float t=t0
        cdef float test
        cdef float pTest
        times=[t0]
        while times[-1]<tMax:
        
            t = t-(1/fmax)*log( float(rand()) / float(RAND_MAX) )
            test = float(rand()) / float(RAND_MAX)
            pTest = PulseIntegrator(t,emin,emax)*areaFrac/fmax
        
            if test <= pTest:
                times.append(t)
        return times


