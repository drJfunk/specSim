from numpy import linspace, arange, empty
from numpy.random import uniform, poisson
from model2counts import model2counts

def PoissonRate(mean,timeBin,dt):

    numPhotons =  poisson(mean,1)
   
    counts = timeBin + uniform(low=0,high=dt,size=numPhotons)
    

    return counts


class photonBuilder(object):

    def __init__(self,tstart,tstop,dt,drm):

        self.chans = arange(0,127,1)
        #self.events =  [[]]*128 #empty chans
        self.events = []
        self.chans = []
        self._CreateTime(tstart,tstop,dt)
        self.drm = drm


    def _EventsFromCounts(self,t):

        
        
        #events = []
        
        for i in arange(0,127,1):

            #i=int(i)
            evts = PoissonRate(self.countSpec[i]*self.dt, t, self.dt)
            #print i
            #print evts
            
            self.events.extend(evts)
            self.chans.extend([i]*len(evts))

        #self.events.append(events)


    def GetModel(self, x):

        return self.model(x, *self.params)

        
            
    def _SetParams(self, params):

        self.params = params


    def SetEvolution(self,paramEvo):
        '''
        This sets the array of paramter evolution used for generating 
        the source counts

        '''

        self.paramEvo = paramEvo

    
    def SetModel(self,model):


        self.model = model


    

    def _CreateTime(self,tstart, tstop, dt):

        self.dt=dt
        self.tstart = tstart
        self.tstop = tstop

        self.timeRange = arange(tstart,tstop,dt)


    def GenerateCounts(self):
        
        for t,p in zip(self.timeRange, self.paramEvo):
            
            self._SetParams(p)
            
            tmpModel = lambda x: self.model(x, *self.params)
            
            m2c = model2counts(self.drm, tmpModel)

            self.countSpec = m2c.counts
            self._EventsFromCounts(t)

            
            


    

