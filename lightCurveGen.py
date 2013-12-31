from bkGroundBuilder import bkGroundBuilder
from photonBuilder import photonBuilder
from rsp import rsp
from tteBuilder import tteBuilder
from numpy import array, argsort, save, load 


class lightCurveGen(object):


    def __init__(self, rspFile, dt, tz, fileName):
        
        self.drm = rsp(rspFile)
        self.source = False
        self.dt =dt
        self.fileName =fileName
        self.tz =tz
        self.bkgFlag=False
        self.sgnFlag=False

    def ConstructBackground(self,tstart,tstop,A,index):



        bkg = bkGroundBuilder(tstart,tstop,self.dt,self.drm)
        bkg.SetSpectrum(A,index)
        bkg.GenerateCounts()

        self.bkgEvents = bkg.events
        self.bkgChans = bkg.chans

        self.t=bkg.timeRange
        self.bkgFlag = True
        del bkg



    def ConstructSignal(self,tstop, model, paramEvo,intType='quad'):
        '''
        The function constructs the source or signal pulse 
        for the simulation. TSTART and DT are taken from the constructor
        but the TSTOP value must be specified here. The param EVO should be
        an array of parameters for EACH parameter in the model. It must
        have as many entries as there are time bins in the sim.

        '''
        sgn = photonBuilder(self.tz,tstop,self.dt,self.drm,intType)
        sgn.SetEvolution(paramEvo)
        
        sgn.SetModel(model)
        sgn.GenerateCounts()
        
        self.sgnEvents = sgn.events
        self.sgnChans = sgn.chans
        self.sgnFlag =True



   
    def SaveBackGround(self,fileName):
        '''
        Load a background. Good for having the same background between 
        files
        '''


        if self.bkgFlag:
            arr= array(zip(self.bkgEvents, self.bkgChans))
            save(fileName,arr)

        else:
            print "No Background created!"


    def LoadBackGround(self,fileName):

        arr = load(fileName)

        self.bkgEvents = arr[:,0]
        self.bkgChans = arr[:,1]

        

    
        
    



    def FormatEvtsChans(self):

        
        self.evts = []
        self.chans = []


        if self.bkgFlag:
            self.evts.extend(self.bkgEvents)
            self.chans.extend(self.bkgChans)
       
       
        if self.sgnFlag:
            self.evts.extend(self.sgnEvents)
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

