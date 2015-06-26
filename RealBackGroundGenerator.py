import astropy.io.fits as fits
import numpy as np
from bisect import bisect
from numpy.random import rand






class RealBackGroundGenerator(object):



    def __init__(self,bakFile):

        bak = fits.open(bakFile)
        rate = bak["SPECTRUM"].data["RATE"]

        self.totalRate = rate.sum()
        self.normRate  = rate/self.totalRate
        self.cumRate   = self.normRate.cumsum()
        self.normedTotal = self.normRate.sum() #Should be one... but check

    
    def generateChannels(self,size):
    
        rv = self.normedTotal*rand(size)
    
        self.chs = np.array(map(lambda x: bisect(self.cumRate,x) , rv)) 
    
        
    def getChannels(self):
        return self.chs

    def GetTotalRate(self):
        return self.totalRate
