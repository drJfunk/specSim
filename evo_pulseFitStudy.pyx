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



#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        
#####################################################################################################        

###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 5.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax





@cython.cdivision(True)
cdef float KRL(float t, float tmax, float c, float r, float d, float fmax):

    cdef float f

    f = fmax*pow((((t+c)/(tmax+c))),r)
    f/=pow(((d+(r*pow((((t+c)/(tmax+c))),(1+r))))/(d+r)),((d+r)/(1+r)))
    
    return f

decay = 1.
#decay = 1.5
#decay = 2.
#decay = 2.5
#decay = 3.
#decay = 3.5
#decay = 4.




@cython.cdivision(True)
cpdef float evo(float ene, float t, p):
    
    
    cdef float A

    maxFlux = 1.
    
    norm = KRL(t,5,0,1.,p[0],maxFlux)
    

    A=p[1]
    Ep=300.
    alpha = -1.
    beta  = -2.2

        
    cdef float val = norm*Band(ene,A,Ep,alpha,beta)
    
    return val




