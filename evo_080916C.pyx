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

        val = A*( pow(x/100., alpha) * exp(-x*(2+alpha)/Ep) )
    else:
        val = A* ( pow( (alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return val



params=np.loadtxt("/Users/jburgess/Research/specSim/grb080916_params.txt")

tStart = params[:,0]+0.1
tStop  = params[:,1]+0.1
norm   = params[:,2]
Ep     = params[:,3]
alpha  = params[:,4]
beta   = params[:,5]

minT = tStart.min()
maxT = tStop.max()



eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t,p):



    cdef float val
    cdef int indx #the index to select from the arrays


    #if t> maxT or t<minT:
    #    return 0.
    
    indx = tStop.searchsorted(t)

    val = 150.*Band(ene,norm[indx],Ep[indx],alpha[indx],beta[indx])

    
        
    return val




