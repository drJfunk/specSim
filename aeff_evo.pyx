
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



cdef float  Band( float x, float A, float Ep, float alpha, float beta):

    cdef float cond = (alpha-beta)*Ep/(2+alpha)
    cdef float val
    if x < cond:

        val = A*( pow(x/100., alpha) * exp(-x*(2+alpha)/Ep) )
    else:
        val = A* ( pow( (alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return val





@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax,params):

        cdef float val
        val = float(quad(evo,emin,emax,args=(t,params))[0])
        return val
        




Eps=np.loadtxt("/Users/jburgess/Research/specSim/Eps_aeffTest.txt")
numIntervals = len(Eps)
intervals = range(numIntervals)






#interpTable = np.loadtxt("/Users/jburgess/Research/specSim/epflux.txt")

#renorm=interp1d(interpTable[:,0],interpTable[:,1])






###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


def pht(float Ep):

    cdef float val
   
    
    b = lambda x: Band(x,1.,Ep,-1.,-2.2) 
    
    val = quad(b,emin,emax)[0]
    
    return val



        
EpGrid = np.logspace(-3.,4.3010299956639813,1000)        
   
flux =[]

for ep in EpGrid:
    flux.append(pht(ep))
    
renorm = interp1d(EpGrid,flux)




@cython.cdivision(True)
cpdef float evo(float ene, float t,p):






    cdef float val
    i= int(floor(t))

    if i > numIntervals-1:
        i=numIntervals-1
    
     
    
    cdef float A


    
    cdef float Ep = Eps[i]
    
    A = 1.E3
    #cdef float A= maxFlux*1.



        
    val = Band(ene,A,Ep,-1.,-2.2)/renorm(Ep)
    


        
    return val




