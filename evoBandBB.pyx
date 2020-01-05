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




eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax




@cython.cdivision(True)
cpdef float  evo(float ene, float t, p):
    
    
    cdef float A   = p[0]
    cdef float Abb = p[1]


        
    cdef float Ep = 300.
    cdef float beta = -2.2
    cdef float alpha = -1.
    cdef float kT = 40.    
        





    cdef float val = Band(ene,A,Ep,alpha,beta)
    
    
    val += BlackBody(ene,Abb,kT)
    
    
        
    return val

