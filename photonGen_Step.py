from scipy.stats import uniform
from scipy.integrate import quad
from numpy import linspace, arange, argsort, array, zeros, ones, mean
from math import log
from photonGen import photonGen
from mnSpecFit.RSPconvolve import RSPconvolve


from numpy import linspace, arange, empty
from numpy.random import uniform, poisson


from evo import evo
import sampler

import cython

def PoissonRate(mean,timeBin,dt):

    numPhotons =  poisson(mean,1)
   
    counts = timeBin + uniform(low=0,high=dt,size=numPhotons)
    

    return counts


class photonGen_Step(photonGen):


    def SetSourceRates(self,rates):
        '''
        the integrated source rates are set here so that it can be
        done faster. Sloppy.... yes


        '''


        self._sourceRates = rates


    def SetRSP(self,rsp):

        self._rsp = rsp
        
    def Set_DT(self,dt):
        '''
        Set the time width for the sim
        '''

        self._dt = dt
    

    def _CreateSourceCurve(self):


        ##This has been modifed to implement cython 2/2/2014
        #self._IntegratePulse()

        #t = arange(self.sourceStart,self.sourceStop,.01)
        #p = map(self._pulse,t)
        #maxFlux = max(p)
        #print "Light curve integrated!"



        #We generate the source enegies first

        
        rspC = RSPconvolve(self._rsp)
        
        


        

        self.eTags=[]
        self.sourceTimes = []
        for t in arange(self.sourceStart,self.sourceStop,self._dt):

            meanT = mean([t,t+self._dt])
            

            
#            tmpModel = lambda e: evo(t,e)


            tmpCounts = zeros(len(rspC.photonE))
            
            #Low res bins

            lowRes = array(map(lambda e: evo(e,meanT,self._additionParams),rspC.lowEne)) 



            medRes = array(map(lambda x:  sum( map(lambda e: evo(e,meanT,self._additionParams),x))/3.,rspC.medEne))

            hiRes =  array(map(lambda x:  sum( map(lambda e:  evo(e,meanT,self._additionParams),x))/7.,rspC.highEne))


            tmpCounts[rspC.lowEval]=lowRes * self._dt
            tmpCounts[rspC.medEval]=medRes * self._dt
            tmpCounts[rspC.highEval]=hiRes * self._dt
            
    
            #print tmpCounts
            rspC.SetModelVec(tmpCounts)
            rspC.CreateModelVector()
            channelRates = array(rspC.GetCounts())[0]
            chNum = 0
            
            for cr in channelRates:
                            
                
                #newPhts = array(self._homoGenSourceCython(t,t+self._dt,cr))#*self.areaFrac)

                # trying to see what the fuck is going on with the rates!!
                #newPhts = PoissonRate(cr*self.areaFrac,t,self._dt)
                newPhts = PoissonRate(cr,t,self._dt) 
                #print newPhts
                #                newPhts = newPhts[newPhts < t+self._dt ]
                
                self.sourceTimes.extend(newPhts)
                self.eTags.extend(ones(len(newPhts))*chNum)

                chNum += 1
                #self.sourceTimes = array(self.sourceTimes)
        print "Generated %d photons from the source"%len(self.sourceTimes)
        
        print "Distributing source photons in energy"
####This is the non parallel version



#        self.srcEnergy = sampler.energize(self.sourceTimes,self.emin,self.emax)


            
        print "Source created!\n\n"



            

