import astropy.io.fits as pf
from numpy import zeros, array, matrix, abs
import math

class rsp(object):

    def __init__(self,rspFileName):


        rspFile = pf.open(rspFileName)
    
        angle = rspFile[2].header["DET_ANG"]
        print "Opening "+rspFileName
        print "Read RSP file with angle: %lf"%angle
        
        

        self.fileName = rspFileName.split('/')[-1]

        


        self.chanData = rspFile[1].data
        self.numEnergyBins = rspFile[2].header['NUMEBINS']
        self.numDetChans = rspFile[2].header['DETCHANS']
        self.det = rspFile[0].header['DETNAM']
        
        self.beta = angle
        self.photonE = (rspFile[2].data["ENERG_LO"], rspFile[2].data["ENERG_HI"])
	self.channelE = (rspFile[1].data["E_MIN"], rspFile[1].data["E_MAX"])
        
        #Main component of object
        self.drm = zeros((self.numEnergyBins, self.numDetChans))

        self._ConstructDRM(rspFile)
        self.binWidths = rspFile[1].data["E_MAX"]-rspFile[1].data["E_MIN"]

        if self.beta:
            self._CalculateProb()


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
                
                #if fc == 0 and nc == self.numDetChans:
                #    self.drm[i]=mData["MATRIX"][i]
                #    print "Uncompressed"
                #else:
                self.drm[i,fc-1:fc+nc]=mData["MATRIX"][i][colIndx:colIndx+nc]
                colIndx+=nc
        self.drm=matrix(self.drm)
        #self.drm=self.drm.transpose()
        del mData

        
    def _CalculateProb(self):
            '''
            Convert the drm into a probablility. To do so, we need the angle 
            relative to the normal of the detector. For the NaI the 
            angle should be relative to the face of the cylinder. For the BGO
            it is given relative to the side of the detector.
            '''

            area = lambda r, h, b: math.fabs(math.pi *(r**2) * math.cos(math.radians(b))) + math.fabs(2 * r * h * math.sin(math.radians(b)))

            if self.det[0].lower() == "n":
                    # NaI's are cylinders of diameter 12.7 cm and thickness 1.27 cm
                    geoArea = area(12.7/2., 1.27, self.beta)
            else:
                    # BGO's are cylinders of diameter 12.7 cm and thickness 12.7 cm.
                    # For the BGO, we have to add an offset of 90 deg to the angle. This
                    # is because for the NaI the angles are calculate relative to the face
                    # of the cylinder whereas for the BGO they are calculated relative to the 
                    # cylindrical face.
                    geoArea = area(12.7/2., 12.7, self.beta + 90.)
            self.geoArea = geoArea
            self.prob = self.drm/geoArea
            # print "Det: %s, Beta: %.1f
    def _FindNearestPhotonBin(self, e,perf=False):
        ''' 
       	Take an energy (e) and find return the 
        corresponding column from the DRM which 
        has been normalised by the geometric area
        '''
        
        idx=(abs(self.photonE[0]-e)).argmin()
        if perf:
            return self.drm[idx,:]
        else:    
            return self.prob[idx,:]
