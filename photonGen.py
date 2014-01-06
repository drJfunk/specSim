from scipy.stats import uniform
from scipy.integrate import quad, quadrature
from numpy import linspace, arange, log

class photonGen(object):

    def __init__(self, bkgStart, bkgStop  ,sourceStart, sourceStop, emin, emax):

        self.srcEnergy = []
        self.bkgEnergy =[]
        self.tStop = bkgStop
        self.tStart = bkgStart
        self.sourceStart = sourceStart
        self.sourceStop = sourceStop
        self.emin = emin
        self.emax = emax


        #self.evts = []



    


    def SetEvolution(self,evo):
        '''
        Pas a fucntion of the form

        f(energy,time) that defines the time and energy 
        evolution of the pulse

        '''
        
        self._specEvo = evo

    def _IntegratePulse(self):

        self._pulse = lambda t: quad(self._specEvo,self.emin,self.emax,args=(t))[0]


    def _nonHomoGen(self, t0, tMax, fmax):
        t=t0
        times=[t0]
        while times[-1]<tMax:
        
            t = t-(1/fmax)*log(uniform.rvs())
        
            if uniform.rvs() <= self._pulse(t)/fmax:
                times.append(t)
        self.sourceTimes = times
        


    def _homoGen(self, t0, tMax, rate):
        t=t0
        times=[t0]
        while times[-1]<tMax:
            times.append(times[-1]-(1/rate)*log(uniform.rvs()))
        
        self.bkgTimes =times



    def _Sampler(self,func,fMax,xMin,xMax):

        flag = True
        while flag:
            energyGuess = uniform.rvs(xMin,xMax-xMin)

            if  uniform.rvs(0,fMax) <= func(energyGuess):
                    flag = False
        
        return energyGuess


    def _CreateSourceCurve(self):


        t = arange(self.sourceStart,self.sourceStop,.001)
        p = map(self._pulse,t)
        maxFlux = max(p)

        self._nonHomoGen(self.sourceStart,self.sourceStop,maxFlux)

        for ti in self.sourceTimes:
            b= lambda en: self._specEvo(en,ti)
            self.srcEnergy.append(self._Sampler(b,b(self.emin),self.emin,self.emax))



    def _CreateBkgCurve(self):

        self._homoGen(self.tStart, self.tStop, self.bkgRate)

        for ti in self.bkgTimes:
            self.bkgEnergy.append(self._Sampler(self._bkgSpec,self._bkgSpec(self.emin),self.emin,self.emax))

            
        
    def _FoldEvents(self):

        pass


