from bkGroundBuilder import bkGroundBuilder
from photonBuilder import photonBuilder
from rsp import rsp
from tteBuilder import tteBuilder
from numpy import array, argsort


class lightCurveGen(object):


    def __init__(self, rspFile, dt, tz, fileName):
        
        self.drm = rsp(rspFile)
        self.source = False
        self.dt =dt
        self.fileName =fileName
        self.tz =tz

    def ConstructBackground(self,tstart,tstop,A,index):



        bkg = bkGroundBuilder(tstart,tstop,self.dt,self.drm)
        bkg.SetSpectrum(A,index)
        bkg.GenerateCounts()

        self.bkgEvents = bkg.events
        self.bkgChans = bkg.chans

        self.t=bkg.timeRange

        del bkg



    def ConstructSignal(self,tstop, model, paramEvo):
        
        sgn = photonBuilder(self.tz,tstop,self.dt,self.drm)
        sgn.SetEvolution(paramEvo)
        
        sgn.SetModel(model)
        sgn.GenerateCounts()
        
        self.sgnEvents = sgn.events
        self.sgnChans = sgn.chans



   

    def FormatEvtsChans(self):

        
        self.evts = []
        self.chans = []

        self.evts.extend(self.bkgEvents)
        self.evts.extend(self.sgnEvents)
        
        self.chans.extend(self.bkgChans)
        self.chans.extend(self.sgnChans)

        self.evts = array(self.evts)
        self.chans = array(self.chans)
        
        sortIndx = argsort(self.evts)

        self.evts=self.evts[sortIndx]
        self.chans=self.chans[sortIndx]


        self.simInfo=[self.tz,self.evts.min(),self.evts.max()]
            
            
            


        

    def MakeTTE(self):

        det = self.drm.det
        tteBuilder(self.fileName,self.evts,self.chans,det,self.simInfo)

