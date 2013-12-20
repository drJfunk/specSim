import astropy.io.fits pf
from numpy import zeros

class rsp(object):

    def __init__(self,rspFileName):


        rspFile = pf.open(rspFileName)
    

        self.chanData = rspFile[1].data
        self.numEnergyBins = rspFile[2].header['NUMEBINS']
        self.numDetChans = rspFile[2].header['DETCHANS']
        
        #Main component of object
        self.drm = zeros((self.numEnergyBins, self.numDetChans))

        self._ConstructDRM(rspFile)

        rspFile.close()

        


    def _ConstructDRM(self,fitFile):
        '''
        Construct the drm from the fits file
        '''
    
        
        mData = fitFile[2].data

        
        self.phtBins = array(zip(mData['ENERG_LO'],mData['ENERG_HI']))

        for fcs,ncs,i in zip(mData["F_CHAN"] , mData["N_CHAN"]  ,range(self.numEnergyBins)):
            colIndx = 0
            for fc,nc in zip(fcs,ncs):

                self.drm[i,fc-1:fc+nc]=mData["MATRIX"][i][colIndx:colIndx+nc]
                colIndx+=nc
        
        del mData

        
