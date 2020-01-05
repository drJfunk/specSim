#Define the Band function
cimport cython

from scipy.integrate import quad, quadrature
from pygsl.sf import synchrotron_1 
import pygsl.errors


from numpy import inf

cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)





#### Synchrotron with pygsl

@cython.cdivision(True)
cpdef float TotalSynchrotron(float x, float A, float eCrit):

    cdef float eta = 9.
    cdef float index = 3.6
    cdef float gammaTh = 3.
    cdef float val

    val = quad(Integrand, 1.,inf, args=(x,A,eCrit,eta,index,gammaTh),epsabs=1.49E-8, epsrel= 1.e-5,maxp1=200 )
    val=val/x
    return val
    
@cython.cdivision(True)
cpdef float Integrand( float gamma, float x , float A, float eCrit, float eta, float index, float gammaTh):
    cdef float val
    
    try:
        val = EDist(A,gamma,eta,gammaTh,index) * synchrotron_1(x/(eCrit*gamma*gamma))[0]
    except pygsl.errors.gsl_Error, err:
        #print err
        val = 0.
    return val


@cython.cdivision(True)
cpdef float EDist(float A, float gamma, float eta, float gammaTh, float index):

        cdef float epsilon
        cdef float val
    
        epsilon = (eta/gammaTh)**(2+index)*exp(-(eta/gammaTh))
    

        if gamma<=eta:

            val = A * (gamma/gammaTh)**2 * exp(-(gamma/gammaTh))

        else:

            val = A * epsilon * (gamma/gammaTh)**(-index)

        return val




####Here is the pulse integrator
#### Uses GSL for fast integration










@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax):

        cdef float val
        val = float(quadrature(evo,emin,emax,args=(t))[0], maxp1=200, maxiter=500)
        return val

        





###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import loadtxt


tmp = loadtxt("sync-intrp-grid.txt")

Ep = tmp[:,0]
flux = tmp[:,1]




    
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


    cdef float indx = -1.
    #cdef float alpha = -1.
    #cdef float beta = -2.2
    
    cdef float A

    

    cdef float Ep = 20* pow(1.+t,indx)
    
    cdef float maxFlux = 10000.
    A = KRL(t,.5,.1, 2.,1.,maxFlux)
    
    cdef float renorm = interpF(Ep) #The interpolated Ep to flux ratio
        
    cdef float val = TotalSynchrotron(ene,A,Ep)/renorm
    
    return val




