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
cdef float  BlackBody( float x, float A, float kT):

   
    cdef float val
    cdef float exponent

    exponent = (x/kT) - 1.

    
    val = A * (x**2) / (exp  (exponent) ) 

    return val



@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax):

        cdef float val
        val = float(quad(evo,emin,emax,args=(t))[0])
        return val



#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        

###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import logspace



def pht(float Ep):

    cdef float val
   
    
    b = lambda x: BlackBody(x,1.,Ep) 
    
    val = quad(b,emin,emax)[0]
    
    return val



        
Ep = logspace(-3.,4.3010299956639813,1000)        
   
flux =[]

for ep in Ep:
    flux.append(pht(ep))
    
interpF = interp1d(Ep,flux)



#### Definition of the KRL Pulse
@cython.cdivision(True)
cdef float KRL(float t, float tmax, float c, float r, float d, float fmax):

    cdef float f

    f = fmax*pow((((t+c)/(tmax+c))),r)
    f/=pow(((d+(r*pow((((t+c)/(tmax+c))),(1+r))))/(d+r)),((d+r)/(1+r)))
    
    return f


@cython.cdivision(True)
cpdef float evo(float ene, float t):
    
    
    cdef float A

    cdef float indx = -1.

    
    cdef float kT = 200.* pow(1.+t,indx)
    
    cdef float maxFlux = 10000.
    A = KRL(t,.5,.1, 2.,1.,maxFlux)
    #cdef float A= maxFlux*1.


    cdef float renorm = interpF(Ep) #The interpolated Ep to flux ratio
        
    cdef float val = BlackBody(ene,A,kT)/renorm
    
    return val




