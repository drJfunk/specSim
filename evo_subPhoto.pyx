#Define the Band function
from subPhotoInterp.subPhotoInterp import subPhotoInterp
cimport cython

from scipy.integrate import quad


cdef extern from "math.h":
    float exp(float)

cdef extern from "math.h":
    float pow(float ,float)

cdef extern from "math.h":
    float log(float)







####Here is the pulse integrator
#### Uses GSL for fast integration

@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax,params):

        cdef float val
        val = float(quad(evo,emin,emax,args=(t,params))[0])
        return val







###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax

from scipy.interpolate import interp1d
from numpy import logspace

sub = subPhotoInterp("/Users/jburgess/Research/mnsf/mnSpecFit/models/TableModel_grid500_res8_Nr1.fits",silent=False)


@cython.cdivision(True)
cpdef float evo(float ene, float t, p):
    
    ##This is a time-independent sim 
   
    cdef float val = sub( ene, 2.,200.,50.,.4 )
    
    return val




