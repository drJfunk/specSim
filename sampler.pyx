cimport cython
from numpy import linspace, array
from evo cimport evo, PulseIntegrator

from scipy.interpolate import interp1d

from libc.stdlib cimport rand, RAND_MAX
cdef extern from "math.h":
    float exp(float)



cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)


@cython.cdivision(True)
cpdef float Sample(float t, float fMax, float xMin, float xMax,params):

        cdef int flag = 1
        cdef float energyGuess
        cdef float test
        cdef float fCheck
        #cdef float randnum
        while flag==1:
           
            energyGuess = xMin +  float(rand()) / float(RAND_MAX/(xMax - xMin))
            test = float(rand()) / float(RAND_MAX/(fMax))
            fCheck = evo(energyGuess,t,params)

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
def nonHomoGen(float t0, float tMax, float emin, float emax, float areaFrac,params):
        '''
        Non-homogeneous poisson process generator
        for a given max rate and time range, this function 
        generates time tags sampled from the energy integrated
        lightcurve.

 
       '''
        print
        print "Integrating pulse to find max"

        
        tmpT=linspace(t0,tMax+1.,1000.)
        tmp = array(map(lambda x: PulseIntegrator(x,emin,emax,params),tmpT))*areaFrac
        tmpIntrp = interp1d(tmpT,tmp,bounds_error=False,fill_value=0.)

        cdef fmax = max(tmp)
        print fmax
        print
        print "Integrated pulse. Now sampling like a mofo"
        print
        cdef float t=t0
        cdef float test
        cdef float pTest
        times=[t0]
        flag =True
        while t<tMax:

            
            
            t = t-(1./fmax)*log( float(rand()) / float(RAND_MAX) )
            test = float(rand()) / float(RAND_MAX)
            pTest = PulseIntegrator(t,emin,emax,params)*areaFrac/fmax
            #pTest = tmpIntrp(t)/fmax


            #print test
            #print pTest
            
            #print t
            
            if test <= pTest:
                times.append(t)
#                print t
        return times

