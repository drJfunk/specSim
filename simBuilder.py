from rsp import rsp
from glob import glob
from numpy import asarray, cumsum, abs, array
from photonGen import photonGen
from numba.random.random import uniform
from tteBuilder import tteBuilder
import sys

class simBuilder(object):
  

    def __init__(self, bn, simDets, ext="",fNameExt ="" ,evo2=None):
       '''
       Builds a simulated burst based off of a certain set of
       RSPs. The entire set of RSPs must be present to determine 
       the flux division between detectors.
    
       The bn is provided to identify the bn of the associated RSPs.
    

       INPUTS:
       
       bn: the burst numbe rof the too be simulated GRB to pull responses
       
       simDets: A list detectors to simulate e.g. simDets = ['n1', 'n6', 'b1']
       
       
       bkgTime: List of background start and stop time [start,stop]
       
       sourceTime: List of source start and stop time [start,stop]
       
       eneSpan:  List of the max and min energies for the spectrum [emin, emax]
       
       bkgParam: List of amplitude and index for the background [A,indx]
    
       evo: Function of time and energy desrcibing the spectral evolution
       

       UPDATE 2/2/2014: Heavily modified photonGen to implement Cython


       '''
     

       self.ext = ext #Used to set folder for holding the files
       self.bn = bn
       self.simDetNames = simDets
       self.evo = None
       self.fNameExt = fNameExt
       self.evo2=evo2

       self.eSpanSet = True


      
           
        
       
    


    def SetBkgParams(self, bkgParam):

        self.A, self.indx = bkgParam
        self.bkgParamSet  = True

    def SetBkgTime(self, bkgTime):

        self.bkgStart, self.bkgStop = bkgTime
        self.bkgTimeSet = True

    def SetSourceTime(self, sourceTime):

        self.sourceStart, self.sourceStop = sourceTime
        self.srcTimeSet = True



    def SetEnergyRange(self, eneSpan):
        
        self.emin, self.emax = eneSpan
        self.eSpanSet = True
    

    def Explode(self):
        
        if self.bkgTimeSet and self.srcTimeSet and self.bkgParamSet and self.eSpanSet:

            assert self.bkgStart < self.bkgStop, "Background times are reversed"
            assert self.sourceStart < self.sourceStop, "Source times are reversed"
            assert self.bkgStart < self.sourceStart, "Source starts before background"
            assert self.bkgStop > self.sourceStop,  "Source ends before background"


            self.simInfo = [self.sourceStart, self.bkgStart, self.bkgStop]

            self._ReadRSP()

        


             #Create the photon generators for each of the desired detectors
            self.simDets =[]
            for i in range(len(self.simDetNames)):
                #pg = photonGen(self.bkgStart,self.bkgStop,self.sourceStart,self.sourceStop,self.emin,self.emax)
                pg = photonGen(self.bkgStart,self.bkgStop,self.sourceStart,self.sourceStop)
                self.simDets.append(pg)

            for frac,sd in zip(self.fractionalAreaNai,self.simDets[1:]):

                sd.SetAreaFraction(frac)

            self._GenerateBackgrounds()


            #Here I need a function that correct reduces the rate of the detectors based on their
            #fraction geometric area. Not sure how to do this
            #FIXED 29/1/2014



            #if self.evo2 != None:
            #    for sd in self.simDets:
            #        sd.SetSecondaryEvolution(self.evo2)


            self._GenerateSignals()


            self._CreateTTE()


        else:
            print "--------> ERROR!!!"
            if not self.bkgTimeSet:
                print "No background times set!"

            if not self.srcTimeSet:
                print "No source times set!"
            
            if not self.bkgParamSet:
                print "Background spectral parameters not set!"

            if not self.eSpanSet:
                print "Max and Min energies not set!"




    def _ReadRSP(self):
        '''
        This method imports the RSPs needed for the simulation
        All fourteen detectors need to be present so that the simulation
        can read calculate the effective area properly


        The effective area will be calculated and then the fractional 
        area will be used to augment the rate for each spectrum


        '''
        
        rsps = glob(self.ext+"*"+self.bn+"*.rsp")
        
        #if len(rsps) != 14:
        #    print "No detectors for "+self.bn
        #    sys.exit()
        #    return

        simRsps = []

        

        self.simDetNames.sort()

        for sd in self.simDetNames:
            try:
                simRsps.append( filter(lambda x: '_'+sd+'_' in x, rsps)[0] )
            except IndexError:
                pass

        

        self.simRsps = map(rsp,simRsps)

        
        bgos = self.simRsps[0]
        nais = self.simRsps[1:]

        

        #totalAreaNai = asarray(map(lambda x: x.geoArea, nais)).sum()
        #print "Total NaI Area: %lf"%totalAreaNai
        
        #totalAreaBgo = asarray(map(lambda x: x.geoArea, bgos)).sum()
        #print "Total BGO Area: %lf"%totalAreaBgo


        self.fractionalAreaNai = map(lambda x: x.geoArea/bgos.geoArea, nais)
        




    def _GenerateBackgrounds(self):
        '''
        For each detector a background will be generated from the master
        background specified

        '''

        for x in self.simDets:

            x.SetBkgParams(self.A,self.indx)
        print "Backgrounda generated"
        


    def _GenerateSignals(self):
        '''
        For each detector generate the pulse from the master parameters
        '''
        for x in self.simDets:

            #x.SetEvolution(self.evo)
            x._CreateSourceCurve()
        print "Signal generated"


    def _CreateTTE(self):

        for x,y,z in zip(self.simDets, self.simRsps, self.simDetNames):

            lc = x.GetLightCurve()
            evts, chans = self._FoldTags(lc,y)

            
            fn = self.ext+z+self.fNameExt+"_simTTE.fit"
            


            tte = tteBuilder(fn, evts, chans, y.det, self.simInfo, y.fileName)
            


    def _FoldTags(self,lc,drm):
        '''
	Loop over input energies and multiply photons
	by the response and obtain a probablility
	distribution in observed energy for that photon.
	'''
	chans = []
	ts = []
	c = 0
	
	tCount, eCount = [], []
	print "Folding time tags through response"
	for tag in lc:
		prob = drm._FindNearestPhotonBin(tag[1])
		pTot = prob.sum()
		
		if pTot == 0:
			continue

		normProb = prob/pTot	

		while pTot > 0:
			r = uniform(0.,1.)
			if r < pTot:
                            r2 = uniform(0.,1.)
                            pCum = cumsum(normProb)
                            idx=(abs(pCum-r2)).argmin() 
                            tCount.append(tag[0])	
                            eCount.append(idx)
			pTot -= 1
		
	
        return array(tCount), array(eCount)
        
