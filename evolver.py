import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


class Evolver(object):
    
    
    def __init__(self):
        
        
        self._alpha = -1.
        self._breakEp = 500.
        self._breakTime = 1.
        
        self._Epindx1 = 0.
        self._Epindx2 = -1.
        
        #KRL
        self._peakTime = 1.
        self._peakNorm = 1.
        self._rise = 1.
        self._decay= 2.
    
    
    def SetPulse(self,peakTime,peakNorm,rise,decay):
        '''
        Set up the pulse profile paramters
        
        - **parameters**, **types**, **return** and **return types**::
        
        :param peakTime: The peak time of the pulse
        :param peakNorm: The Normalization of the pulse at peak time
        :param rise: The power law rise index
        :param decay: The power law decay index
        '''
        self._peakTime = float(peakTime)
        self._peakNorm = float(peakNorm)
        self._rise     = float(rise)
        self._decay    = float(decay)
        
    def SetSpectrum(self,Ep,alpha):
        '''
        Set up the spectral parameters
        
        - **parameters**, **types**, **return** and **return types**::
        
        :param Ep: Ep and the time of the evolution break
        :param alpha: The low energy spectral index
        '''
        
        self._breakEp = float(Ep)
        self._alpha   = float(alpha)
    
    def SetSpectralEvolution(self,breakTime,EpIndex1,EpIndex2):
        '''
        Set up the spectral parameters
        
        - **parameters**, **types**, **return** and **return types**::
        
        :param breakTime: When Ep breaks
        :param EpIndex1: Ep index before the break
        :param EpIndex2: Ep index after the break
        '''
        
        self._breakTime = float(breakTime)
        self._Epindx1  = float(EpIndex1)
        self._Epindx2  = float(EpIndex2)
            
        
    
    
    
    def IntegratedSpectrum(self,energy,t0,tf):
        '''
        Returns the specific flux at energy integrated betweeen t0 and tf
        
        - **parameters**, **types**, **return** and **return types**::
        
        :param energy: evaluation energy
        :param t0: start of the temporal integration
        :param t0: end of the temporal integration
        :returns: temporally integrated specfic photon flux
        '''
        
        
        return quad(self._cutoffPL,t0,tf,args=(energy) ,epsabs=1.49e-19 )[0]
    
    
    
    def PlotPulse(self,t0,tf):
        fig = plt.figure(100)
        pulseax = fig.add_subplot(111)
        epax    = pulseax.twinx()
        
        
        timeGrid = np.linspace(t0,tf,1000)
        
        pulseax.plot(timeGrid,self._KRL(timeGrid,self._peakTime,0.,self._rise,self._decay,self._peakNorm),'b')
        
        epax.semilogy(timeGrid,self._brokenPL(timeGrid),'r')
        
        pulseax.set_xlabel('Time')
        pulseax.set_ylabel('Flux')
        epax.set_ylabel('Ep')
        
        
    def PlotVFv(self,t0,tf,emin=10.,emax=1E4):
        
        
        fig = plt.figure(200)
        ax = fig.add_subplot(111)
        
        
        energyGrid = np.logspace(np.log10(emin),np.log10(emax),100)
        
        ax.loglog(energyGrid,energyGrid**2 * np.array([self.IntegratedSpectrum(e,t0,tf) for e in energyGrid]))
        
        
        
        ax.set_xlabel('Energy')
        ax.set_ylabel('vFv')
        ax.set_ylim(bottom=1E2)
    
    
    def _brokenPL(self,time):
        time = np.array(time,ndmin=1,copy=False)
       
        vals = np.zeros(time.flatten().shape[0])
        
        cond = time<self._breakTime
        #print cond
        vals[cond] = (time[cond]/self._breakTime)**self._Epindx1
        vals[~cond] = (time[~cond]/self._breakTime)**self._Epindx2
        return self._breakEp*vals
        

    def _cutoffPL(self,time,energy):
    
        N=self._KRL(time,self._peakTime,0.,self._rise,self._decay,self._peakNorm)
            
        return N*(energy/100.)**self._alpha * np.exp(-energy/self._brokenPL(time))

    def _KRL(self,t,tmax,c,r,d,fmax):

        f = (fmax*(np.power((((t+c)/(tmax+c))),r)/np.power(((d+(r*np.power((((t+c)/(tmax+c))),(1+r))))/(d+r)),((d+r)/(1+r)))))
        return f    
