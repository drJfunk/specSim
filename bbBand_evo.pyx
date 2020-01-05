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


@cython.cdivision(True)
cdef float  Band( float x, float A, float Ep, float alpha, float beta):

    cdef float cond = (alpha-beta)*Ep/(2+alpha)
    cdef float val
    if x < cond:

        val = A*( pow(x/100., alpha) * exp(-x*(2+alpha)/Ep) )
    else:
        val = A* ( pow( (alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*pow(x/100.,beta))

    return val




eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import logspace

cdef float alpha = -1.
cdef float beta = -2.2


def pht(float Ep):

    cdef float val
   
    
    b = lambda x: Band(x,1.,Ep,alpha,beta) 
    
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

def bbPht(kT):
    
    b = lambda x: BlackBody(x,1.,kT) 
    
    cdef float val = quad(b,emin,emax)[0]
    
    return val
kT = logspace(-4,4,1000)
flux =[]

for kt in kT:
    flux.append(pht(kt))
    
interpBB = interp1d(kT,flux)

@cython.cdivision(True)
cpdef float  evo(float ene, float t):
    
    
    cdef float A
    
    cdef float indx = -1.

    
    cdef float Ep = 2000.* pow(1.+t,indx)
    
    cdef float maxFlux = 10000.
    A = KRL(t,.5,.1, 2.,1.,maxFlux)
   


    cdef float renorm = interpF(Ep) #The interpolated Ep to flux ratio
        
    cdef float val = Band(ene,A,Ep,alpha,beta)/renorm
    
    
    ######BB
    
    cdef float tbreak = 5.
    
    cdef float indx1 = -.66
    cdef float indx2 = -.7
    cdef float kT
    if t<= tbreak:
        kT = 80.*pow(1+t,indx1)
    else:
        kT= 80*pow(tbreak,indx1-indx2)*pow(1+t,indx2)
        
  
           
    maxFlux = 1E-0
    A = KRL(t,.5,.1, 2.,1.,maxFlux)
    renorm = interpBB(kT)
    val += BlackBody(ene,A,kT)/renorm
    
    
        
    return val

