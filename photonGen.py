from scipy.stats import uniform
from scipy.integrate import quad
from numpy import linspace, arange, argsort, array
from math import log
from RealBackGroundGenerator import RealBackGroundGenerator

#from numba.random.random import uniform
#from Tkinter import *
#from progressBar import Meter

import evo
import sampler

import cython



class photonGen(object):

    def __init__(self, bkgStart, bkgStop  ,sourceStart, sourceStop, emin=None, emax=None,noBkg=False):
        '''
        This class creates an object that performs the task of 
        generating photons time for both a pulse and a background.
        The photons are distributed in energy according to a time-dependent
        coded spectrum.

        INPUTS:
        bkgStart: Starting time of the background
        bkgStop: Stopign time of the background

        sourceStart: Starting time of the sourceStart
        sourceStop: Stoping time of the source

        emin: The minimum energy of the input spectrum
        emax: The maximum energy of the input spectrum


        The call structure is as follows:

        Instantiate the class:
          pGen = photonGen(**args)


        Input the background spectral params:
          pGen.SetBkgParams(**args)


        A function must be created that defines the spectral evolution.
        As of now, the function must be analytic but plans for reading a
        txt file will be implemented. The function for evo is coded and input

           pGen.SetEvolution(evo)

        The object now computes the energy integrated lightcurves
        for both the source and background and then samples these 
        lightcurves for time tags. These time tags then go back in
        so that the time dependent spectrum can be sampled.

        The result will be a list of time tags with photon energies
        that will be fed to the response code for calculating counts

        UPDATE 2/2/2014: Heavily modified to implemenet Cython bindings
        ADD MORE DETAIL


        '''

        self.realBkg = False
        self.noBkg = noBkg
        self.srcEnergy = []
        self.bkgEnergy =[]
        self.tStop = bkgStop
        self.tStart = bkgStart
        self.sourceStart = sourceStart
        self.sourceStop = sourceStop
        self.emin = evo.eMin
        self.emax = evo.eMax
        self.areaFrac = 1.

        self._evo2 = False


        
    


    def SetAreaFraction(self, frac):
        '''
        Sets the fraction of the effective area to a BGO
        detector so that the rate is modified.

        '''

        self.areaFrac = frac
        print self.areaFrac


    def SetBkgParams(self,amp,index):
        '''
        The background is modeled as 
        power law with input amplitude and
        index. This method sets those parameters
        and then calls method to generate the background 
        lightcurve
        

        '''



        self.bkgAmp = amp
        self.bkgIndex = index
        

        self._CreateBkgCurve()

    def SetSecondaryEvolution(self,evo):
        '''
        In the case that the evolution function is weighted or hard to 
        calculate, this allows for a secondary evolution to be set that
        simpler to compute. 


        '''
        self._evo2 = evo



    def SetEvolution(self,evo):
        '''
        Pass a fucntion of the form

        f(energy,time) that defines the time and energy 
        evolution of the pulse

        The method then calls the fucntion to create the 
        source lightcurve



        '''
        
        self._specEvo = evo

        self._CreateSourceCurve()

            


    def _IntegratePulse(self):
        '''
        In order to have the correct number of photons in the lightcurve
        the spectral evolution must be integrated of energy to get the 
        total possion rate. This method implemenets that process.

        '''


        self._pulse = lambda t: quad(self._specEvo,self.emin,self.emax,args=(t,self._additionParams))[0]*self.areaFrac
        



    def _nonHomoGenCython(self, t0, tMax):

        print "Invoking Cython Non-Homogeneus Sampler"
       
        self.sourceTimes = sampler.nonHomoGen(t0,tMax,self.emin,self.emax,self.areaFrac,self._additionParams)



    def _nonHomoGen(self, t0, tMax, fmax):
        '''
        Non-homogeneous poisson process generator
        for a given max rate and time range, this function 
        generates time tags sampled from the energy integrated
        lightcurve.

        '''


        print "Invoking python version of non-homo gen"

        t=t0
        times=[t0]
        while times[-1]<tMax:
        
            t = t-(1/fmax)*log(uniform.rvs(0.,1.))
        
            if uniform.rvs(0.,1.) <= self._pulse(t)/fmax:
                times.append(t)
        self.sourceTimes = times
        print "There were %d photons generated.\nDistributing in energy"%len(self.sourceTimes)
        


    def _homoGenCython(self, t0, tMax, rate):


        self.bkgTimes = sampler.homoGen(t0,tMax,rate)

    def _homoGenSourceCython(self, t0, tMax, rate):


        return sampler.homoGen(t0,tMax,rate)


    def _homoGen(self, t0, tMax, rate):
        t=t0
        times=[t0]
        while times[-1]<tMax:
            times.append(times[-1]-(1/rate)*log(uniform.rvs(0.,1.)))
        
        self.bkgTimes =times


    def _SamplerInv(self,xMin,xMax):
        '''
        Since the bkg is a power-law, we can use an inverse transform 
        sampler which is faster.

        '''

        indx = self.bkgIndex

        #pick a random number
        r =  uniform.rvs(0.,1.)

        #For a power law, the inverse CDF is easy to calculate
        energy = ((xMax**(indx+1)-xMin**(indx+1))*r + xMin**(indx+1))**(1./(indx+1))

        return energy
        


    def _SetParams(self, params):

        self._additionParams = params

    def _SamplerRej(self,func,fMax,xMin,xMax):

        flag = True
        while flag:
            energyGuess = uniform.rvs(xMin,xMax-xMin)

            if  uniform.rvs(0,fMax) <= func(energyGuess):
                    flag = False
        
        return energyGuess


    def _CreateSourceCurve(self):


        ##This has been modifed to implement cython 2/2/2014
        #self._IntegratePulse()

        #t = arange(self.sourceStart,self.sourceStop,.01)
        #p = map(self._pulse,t)
        #maxFlux = max(p)
        #print "Light curve integrated!"

        self._nonHomoGenCython(self.sourceStart,self.sourceStop)

        print "Generated %d photons from the source"%len(self.sourceTimes)
        


        print "Distributing source photons in energy"
####This is the non parallel version
        for ti,i in zip(self.sourceTimes,xrange(len(self.sourceTimes))):

            fMax = evo.evo(self.emin,ti,self._additionParams)
            self.srcEnergy.append(sampler.Sample(ti,fMax,self.emin,self.emax,self._additionParams))


            
        print "Source created!\n\n"


    def _bkgSpec(self,ene):

        return self.bkgAmp*(ene)**(self.bkgIndex)


    def _IntegrateBkg(self):


        if self.realBkg:
            print
            print "Getting bkg rate from BAK file"
            self.bkgRate = self.bkgGen.GetTotalRate()
        else:
            self.bkgRate = quad(self._bkgSpec, self.emin,self.emax)[0]*self.areaFrac
        #print self.bkgRate
        #print self.areaFrac

    def SetBakFile(self,bakFile):
        self.realBkg = True

        self.bkgGen = RealBackGroundGenerator(bakFile)
        
    def _CreateBkgCurve(self):

        if self.noBkg:
            self.bkgTimes = []
            self.bkgEnergy = []
            print "No Background Created"
            return
        
        
        self._IntegrateBkg()

        self._homoGen(self.tStart, self.tStop, self.bkgRate)
        
        print "Generated %d photons from the background"%len(self.bkgTimes)
        print "Distributing background photons in energy"

        if self.realBkg:
            self.bkgGen.generateChannels(len(self.bkgTimes))

        else:
            for ti in self.bkgTimes:
                self.bkgEnergy.append(self._SamplerInv(self.emin,self.emax))
        print "Background created"

            
    def GetLightCurve(self):
        '''
        Return the time sorted combined source and 
        background photon time tags

        '''

        # combine the source and background times
        tags = []
        tags.extend(self.bkgTimes)
        tags.extend(self.sourceTimes)

        tags = array(tags)
        ene = []

        # combine the source and background energies
        ene.extend(self.bkgEnergy)
        ene.extend(self.srcEnergy)

        ene = array(ene)

        #Sort the times and return an array of sort inidices 
        #to sort the energies
        indx = argsort(tags)

        tags = tags[indx]
        ene = ene[indx]

        lc = array(zip(tags,ene))

        
        
        return lc
