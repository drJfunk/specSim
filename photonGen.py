from scipy.stats import uniform
from scipy.integrate import quad, quadrature


class photonGen(object):

    def __init__(self):

        self.srcEnergy = []
        self.bkgEnergy =[]
        #self.evts = []



    def _specEvo(self, t, ene):
        pass


    def _pulse(self, t):
        pass

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



    def _Sampler(func,fMax,xMin,xMax):

        flag = True
        while flag:
            energyGuess = uniform.rvs(xMin,xMax-xMin)

            if  uniform.rvs(0,fMax) <= func(energyGuess):
                    flag = False
        return energyGuess


    def _CreateSourceCurve(self):


        t = linspace(self.sourceStart,self.sourceStop,.01)
        p = map(self._pulse,t)
        maxFlux = max(p)

        self._nonHomoGen(self.sourceStart,self.sourceStop,maxFlux)

        for ti in self.sourceTimes:
            b= lambda en: self._specEvo(ti,en)
            self.srcEnergy.append(self._Sampler(b,b(self.emin),self.emin,self.emax))



    def _CreateBkgCurve(self):

        self._homoGen(self.tStart, self.tStop, self.bkgRate)

        for ti in self.bkgTimes:
            self.bkgEnergy.append(self._Sampler(self._bkgSpec,self._bkgSpec(self.emin),self.emin,self.emax))

            
        
    def _FoldEvents(self):

        pass


