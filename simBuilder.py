from rsp import rsp
from glob import glob
from numpy import asarray, cumsum, abs, array
from photonGen import photonGen
from numba.random.random import uniform
from tteBuilder import tteBuilder

class simBuilder(object):
  

    def __init__(self, bn, simDets, bkgTime, sourceTime, eneSpan, bkgParam, evo, ext=""):
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
       

       '''
     

       self.ext = ext #Used to set folder for holding the files

       emin, emax = eneSpan

       bkgStart, bkgStop = bkgTime
       
       sourceStart, sourceStop = sourceTime
       

       assert bkgStart < bkgStop, "Background times are reversed"
       assert sourceStart < sourceStop, "Source times are reversed"
       assert bkgStart < sourceStart, "Source starts before background"
       assert bkgStop > sourceStop,  "Source ends before background"
        
       self.A, self.indx = bkgParam
       
       self.bn = bn
       
       self.simDetNames = simDets

       self.simInfo = [sourceStart, bkgStart, bkgStop]

       self.evo = evo

        
       self._ReadRSP()

        


       #Create the photon generators for each of the desired detectors
       self.simDets = [photonGen(bkgStart,bkgStop,sourceStart,sourceStop,emin,emax)] * len(simDets)
      

       self._GenerateBackgrounds()

       
       #Here I need a function that correct reduces the rate of the detectors based on their
       #fraction geometric area. Not sure how to do this
       
       self._GenerateSignals()


       self._CreateTTE()
    


        
    




    def _ReadRSP(self):
        '''
        This method imports the RSPs needed for the simulation
        All fourteen detectors need to be present so that the simulation
        can read calculate the effective area properly


        The effective area will be calculated and then the fractional 
        area will be used to augment the rate for each spectrum


        '''
        
        rsps = glob(self.ext+"*"+self.bn+"*.rsp")
        

        simRsps = []

        for sd in self.simDetNames:
            try:
                simRsps.append( filter(lambda x: '_'+sd+'_' in x, rsps)[0] )
            except IndexError:
                pass

        rsps = map(rsp,rsps)

        self.simRsps = map(rsp,simRsps)

        
        bgos = rsps[:2]
        nais = rsps[2:]

        

        totalAreaNai = asarray(map(lambda x: x.geoArea, nais)).sum()
        print "Total NaI Area: %lf"%totalAreaNai
        
        totalAreaBgo = asarray(map(lambda x: x.geoArea, bgos)).sum()
        print "Total BGO Area: %lf"%totalAreaBgo


        fractionalAreaNai = map(lambda x: x.geoArea/totalAreaNai, nais)
        




    def _GenerateBackgrounds(self):
        '''
        For each detector a background will be generated from the master
        background specified

        '''

        for x in self.simDets:

            x.SetBkgParams(self.A,self.indx)
        print "Background generated"
        


    def _GenerateSignals(self):
        '''
        For each detector generate the pulse from the master parameters
        '''
        for x in self.simDets:

            x.SetEvolution(self.evo)

        print "Signal generated"


    def _CreateTTE(self):

        for x,y,z in zip(self.simDets, self.simRsps, self.simDetNames):

            lc = x.GetLightCurve()
            evts, chans = self._FoldTags(lc,y)

            
            fn = self.ext+z+"_simTTE.fit"
            


            tte = tteBuilder(fn, evts, chans, y.det, self.simInfo)
            


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
				pCum = cumsum(normProb)
				idx=(abs(pCum-r)).argmin() 
				tCount.append(tag[0])	
				eCount.append(idx)
			pTot -= 1
		
	
        return array(tCount), array(eCount)
        
