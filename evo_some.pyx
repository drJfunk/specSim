cimport cython
from math import floor
import numpy as np
#cimport numpy as np
from scipy.interpolate import interp1d
from numpy import *

from pygsl.sf import synchrotron_1
import pygsl.errors

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

        



#for i in range(1000):

def makeGamma(norm,A,gMax,p,DT):
    
    #gMin = 1000.
    N=500
    gMin=1000.
    def Source(x):
    
        if (x<gMin or x> gMax):
            return 0.
        else:
            return power(x,-p)*(1-p)* 1./ (power(gMax,1-p)-pow(gMin,1-p))
    
    step = exp(1./N*log(gMax*1.1))

    gamma  = zeros(N)
    fgamma = zeros(N)
    fgammatp1 =  zeros(N)
    G = zeros(N+1)
    
    factor = 25352.525
    
    cool = A**2*factor
    
    
    
    for i in range(N):

        gamma[i]=power(step,i)
        fgamma[i]=Source(gamma[i])

        if i<N-1:
            G[i] = 0.5*(gamma[i]+gamma[i]*step)
        else:
            G[N-1] = 0.5*(gamma[i]+gamma[i]*step)

    deltaGamma = G[N-1]-G[N-2]
    fgammatp1[N-1]=fgamma[N-1]/(1. + (DT*cool*gamma[N-1]*gamma[N-1] )/ deltaGamma)

    for j in range(N-2,0,-1):
        deltaGamma = 0.5*(G[j]-G[j-1])
        gdotp=cool*gamma[j+1]*gamma[j+1]
        gdotm=cool*gamma[j]*gamma[j]

        V3 = (DT*gdotp )/ deltaGamma
        V2 = 1. + (DT*gdotm )/ deltaGamma
        fgammatp1[j] = (fgamma[j]+Source(gamma[j])*DT + V3*fgammatp1[j+1])/V2

    deltaGamma = G[1]-G[0]
    fgammatp1[0]=fgamma[0]+(DT*cool*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma
    return [norm*(fgammatp1),gamma]

def sync(ene,fgamma,gamma,A):
    
    N=500
    y=zeros(N)
    for i in range(N):
        try:
            y[i] = fgamma[i]*synchrotron_1(ene/(A*gamma[i]*gamma[i]))[0]
        except pygsl.errors.gsl_Error, err:
            pass
    #return y
    s = y[0] + 2.0 * sum(y[1:N-1]) + y[N-1]
    h = float(gamma[-1] - gamma[0]) / N
    #print s
    return A*s * h / (2.0*ene)




eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t, p):



    cdef float val
    
    

    tcool = p[4]/(p[1]*p[1]*25352.525)
    fgamma,gamma = makeGamma(p[0],p[1],p[2],p[3],tcool)

    val = p[5]*sync(ene,fgamma,gamma,p[5]*p[1])
    
    
        
    return val




