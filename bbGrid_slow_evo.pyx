
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


cdef extern from "synchrotron.h":
    double synchrotron(double,double,double,double)
    double synchrotronComplex(double,double,double,double,double,double)
    double synchrotronPL(double, double, double, double, double)
    double synchrotronFast(double, double, double, double, double)
    double synchrotronPL_cutoff(double, double, double, double, double, double)
    double synchrotron_cutoff(double, double, double, double, double)
    double SSC(double, double, double , double )






@cython.cdivision(True)
cpdef float PulseIntegrator(float t, float emin, float emax):
    cdef float val
    val = float(quad(evo,emin,emax,args=(t))[0])
    return val

@cython.cdivision(True)
cpdef float bb(float ene, float A, float kT):

    cdef float val
    val =  A*ene**2.*( exp(ene/kT) -1.)**-1.
    return val



bbGrid = np.loadtxt("/Users/jburgess/Research/specSim/bbGrid_slow_Ep300keV.txt") 
#bbGrid = np.loadtxt("/Users/jburgess/Research/specSim/bbGrid_slow_Ep1MeV.txt")
#bbGrid = np.loadtxt("/Users/jburgess/Research/specSim/bbGrid_TEST.txt")
kT = bbGrid[:,0]
frac = bbGrid[:,1]


numIntervals = len(kT)
intervals = range(numIntervals)

print numIntervals


interpTable = np.loadtxt("/Users/jburgess/Research/specSim/epflux.txt")

renorm=interp1d(interpTable[:,0],interpTable[:,1])






###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax







@cython.cdivision(True)
cpdef float evo(float ene, float t):



    cdef double val
    cdef double synch

    cdef double kt

    i= int(floor(t))

    if i > numIntervals-1:
        i=numIntervals-1

        
    fr=frac[i]
    kt = kT[i]
    
    
    Ep = 300.
    Estar = Ep/(3.*90.**2.)
    #SN = 30.
    #bkg = 265.64
    synch = 100.*synchrotronPL(ene,1.2499E-0,Estar,3.5,90.)/renorm(Estar) 
    
    
    
    val =bb(ene,fr,kt)



    val +=synch
    #val = 1E-2

    return val




