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

        




Eps=np.loadtxt("/Users/jburgess/Research/specSim/Eps.txt")
numIntervals = len(Eps)
intervals = range(numIntervals)

print numIntervals


interpTable = np.loadtxt("/Users/jburgess/Research/specSim/epflux.txt")

renorm=interp1d(interpTable[:,0],interpTable[:,1])





print numIntervals
###### THIS WILL CHANGE DEPENDING ON EVOLUTION DESIRED
eMin = 6.
eMax = 50000.

cdef float emin = eMin
cdef float emax = eMax


@cython.cdivision(True)
cpdef float evo(float ene, float t):



    cdef float val
    i= int(floor(t))

    if i > numIntervals-1:
        i=numIntervals-1
    
    Estar = Eps[i]/(3.*90.**2.)
    SN = 30.
    bkg = 186.73665961
    val = (SN*bkg)*synchrotronPL(ene,1.,Estar,3.5,90.)/renorm(Estar)
        
    return val




